#!/usr/bin/env python3
"""Run isolated CLI correctness checks in the regression runner's temporary cwd."""
import os
import re
from pathlib import Path
import subprocess
import sys

exe = sys.argv[1]

def write(name, content):
    Path(name).write_text(content)

def run(args, success=True):
    result = subprocess.run([exe, *args], capture_output=True, text=True, timeout=10)
    assert result.returncode >= 0, (args, result.returncode, result.stderr)
    assert 'AddressSanitizer' not in result.stderr and 'runtime error:' not in result.stderr, result.stderr
    assert (result.returncode == 0) == success, (args, result.returncode, result.stderr)
    return result

def rows(result):
    lines = result.stdout.splitlines()
    assert lines and lines[0] == 'ID EXPT AVDIFF', result.stdout
    return {fields[0]: float(fields[2]) for fields in (line.split() for line in lines[1:])}

write('exp', 'ID measured\na 0\nb 0\nc 0\n')
write('pred1', 'ID predicted\na 2\nb 1\nc 10\n')
write('pred2', 'ID predicted\na -2\nb 1\nc 10\n')
assert rows(run(['-E', 'exp', '-d', 'pred1', 'pred2'])) == {'a': 2, 'b': 1, 'c': 10}
for count in [1, 2, 3, 4, 100]:
    for balanced in [False, True]:
        args = ['-E', 'exp', '-d', '-w', str(count)]
        if balanced:
            args += ['-b']
        result = rows(run([*args, 'pred1']))
        assert len(result) == min(count, 3), (args, result)
        assert 'c' in result, (args, result)
for predictions in ['ID predicted\na 2\nb 1\nc 10\n',
                    'ID predicted\na -2\nb -1\nc -10\n',
                    'ID predicted\na 2\nb 0\nc -10\n']:
    write('threshold_pred', predictions)
    for threshold, expected in [('0.5', {k for k, v in rows(run(['-E', 'exp', '-d', 'threshold_pred'])).items() if abs(v) >= .5}),
                                ('10', {'c'}), ('20', set())]:
        assert set(rows(run(['-E', 'exp', '-d', '-W', threshold, 'threshold_pred']))) == expected
write('one_exp', 'ID measured\na 0\n')
write('one_pred', 'ID predicted\na 2\n')
assert rows(run(['-E', 'one_exp', '-d', '-w', '1', '-b', 'one_pred'])) == {'a': -2}
write('mixed', 'ID predicted\na -10\nb -9\nc 8\n')
assert set(rows(run(['-E', 'exp', '-d', '-w', '2', '-b', 'mixed']))) == {'a', 'c'}
assert rows(run(['-E', 'exp', '-d', '-R', '^pred[12]$'])) == {'a': 2, 'b': 1, 'c': 10}
run(['-E', 'exp', '-R', '['], success=False)
run(['-E', 'exp', '-R', '^no_matching_files$'], success=False)
write('bad_pred', 'ID predicted\na 2\nb nonsense\nc 10\n')
run(['-E', 'exp', 'bad_pred'], success=False)
assert rows(run(['-E', 'exp', '-d', '-F', 'bad_pred'])) == {'a': 2, 'c': 10}
write('short_pred', 'ID predicted\na 2\nb\nc 10\n')
assert rows(run(['-E', 'exp', '-d', '-F', 'short_pred'])) == {'a': 2, 'c': 10}
write('empty_pred', 'ID predicted\n')
run(['-E', 'exp', 'empty_pred'], success=False)
write('duplicates', 'ID measured\na 0\na 100\n')
run(['-E', 'duplicates', 'one_pred'], success=False)
assert rows(run(['-E', 'duplicates', '-d', '-h', 'one_pred'])) == {'a': 2}
write('class_exp', 'ID class\na A\na B\n')
write('class_pred', 'ID class\na A\n')
assert 'a A 1 1 1' in run(['-C', '-d', '-h', '-E', 'class_exp', 'class_pred']).stdout
write('class_bad', 'ID class\na B\nb\na A\n')
assert 'a A 2 1 0.5' in run(['-C', '-d', '-h', '-F', '-E', 'class_exp', 'class_bad']).stdout
# Reject aliases before opening output; input bytes must remain intact.
original = Path('exp').read_bytes()
for output in ['exp', './exp']:
    run(['-E', 'exp', '-Y', output, 'pred1'], success=False)
    assert Path('exp').read_bytes() == original
os.link('exp', 'hardlink')
run(['-E', 'exp', '-Y', 'hardlink', 'pred1'], success=False)
assert Path('exp').read_bytes() == original
write('prediction.textproto', Path('pred1').read_text())
run(['-E', 'exp', '-T', 'prediction', 'prediction.textproto'], success=False)
assert Path('prediction.textproto').read_text() == Path('pred1').read_text()
write('files.list', 'pred1\npred2\n')
run(['-E', 'exp', '-Y', 'pred1', 'F:files.list'], success=False)
run(['-E', 'exp', '-Y', 'pred1', '-R', '^pred[12]$'], success=False)
run(['-E', 'exp', '-Y', 'same.textproto', '-T', 'same', 'pred1'], success=False)
assert not Path('same.textproto').exists()
# Reversing display must not change either filter's membership, including ties.
for selection in [[], ['-w', '2'], ['-w', '2', '-b'], ['-W', '1.5'], ['-W', '20']]:
    args = ['-E', 'exp', '-d', *selection, 'pred1', 'pred2']
    normal = run(args).stdout.splitlines()
    reverse = run(['-r', *args]).stdout.splitlines()
    assert reverse == normal[:1] + normal[:0:-1], (args, normal, reverse)
write('tie_exp', 'ID measured\nz 0\naa 0\na 0\n')
write('tie_pred', 'ID predicted\naa 1\nz 1\na 1\n')
normal = run(['-E', 'tie_exp', '-d', 'tie_pred']).stdout.splitlines()
assert [row.split()[0] for row in normal[1:]] == ['a', 'aa', 'z']
reverse = run(['-E', 'tie_exp', '-d', '-r', 'tie_pred']).stdout.splitlines()
assert reverse == normal[:1] + normal[:0:-1]
selected = run(['-E', 'tie_exp', '-d', '-w', '2', 'tie_pred']).stdout
write('tie_exp', 'ID measured\na 0\nz 0\naa 0\n')
write('tie_pred', 'ID predicted\nz 1\na 1\naa 1\n')
assert run(['-E', 'tie_exp', '-d', 'tie_pred']).stdout.splitlines() == normal
assert run(['-E', 'tie_exp', '-d', '-w', '2', 'tie_pred']).stdout == selected
# Summary and textproto ordering follow the same reversal as primary output.
for suffix, reverse_flag in [('normal', []), ('reverse', ['-r'])]:
    run(['-E', 'exp', '-d', '-Y', f'summary_{suffix}', '-T', f'proto_{suffix}',
         *reverse_flag, 'pred1', 'pred2'])
summary = Path('summary_normal').read_text().splitlines()
assert Path('summary_reverse').read_text().splitlines() == summary[:1] + summary[:0:-1]
assert Path('proto_reverse.textproto').read_text().splitlines() == Path('proto_normal.textproto').read_text().splitlines()[::-1]
write('class_tie_exp', 'ID class\nz A\na A\n')
write('class_tie_pred', 'ID class\na A\nz A\n')
normal = run(['-C', '-E', 'class_tie_exp', '-d', 'class_tie_pred']).stdout.splitlines()
assert [row.split()[0] for row in normal[1:]] == ['a', 'z']
assert run(['-C', '-E', 'class_tie_exp', '-d', '-r', 'class_tie_pred']).stdout.splitlines() == normal[:1] + normal[:0:-1]
# Opposite-sign errors must not cancel in default worst/threshold selection.
assert rows(run(['-E', 'exp', '-d', '-w', '2', 'pred1', 'pred2'])) == {'a': 2, 'c': 10}
assert rows(run(['-E', 'exp', '-d', '-W', '1.5', 'pred1', 'pred2'])) == {'a': 2, 'c': 10}
assert rows(run(['-E', 'exp', '-d', '-s', '-w', '2', 'pred1', 'pred2'])) == {'b': -1, 'c': -10}
assert rows(run(['-E', 'exp', '-d', '-s', '-W', '1.5', 'pred1', 'pred2'])) == {'c': -10}
run(['-E', 'exp', '-b', 'pred1'], success=False)
# Single-observation summaries must avoid variance(), even with assertions enabled.
for flags in [[], ['-s']]:
    run(['-E', 'one_exp', '-Y', 'single_summary', '-T', 'single_proto', *flags, 'one_pred'])
    record = Path('single_summary').read_text().splitlines()[1].split()
    assert record[0] == 'a' and record[2] == '1'
    assert float(record[-1]) == 0.0
    proto = Path('single_proto.textproto').read_text()
    assert 'npred: 1' in proto and 'std_difference: 0' in proto
# Multiple observations retain sample standard deviation of the selected metric.
write('spread1', 'ID predicted\na 2\n')
write('spread2', 'ID predicted\na -4\n')
for flags, mean, stddev in [([], 3.0, 2.0 ** .5), (['-s'], 1.0, 18.0 ** .5)]:
    run(['-E', 'one_exp', '-Y', 'spread_summary', '-T', 'spread_proto', *flags, 'spread1', 'spread2'])
    record = Path('spread_summary').read_text().splitlines()[1].split()
    assert abs(float(record[6]) - mean) < 1e-5
    assert abs(float(record[-1]) - stddev) < 1e-5
    proto = Path('spread_proto.textproto').read_text()
    value = float(re.search(r'std_difference: ([^ ]+)', proto).group(1))
    assert abs(value - stddev) < 1e-5
# Apply minimum valid-record count before selection and across all regression outputs.
write('coverage1', 'ID predicted\na 2\nb 1\nc 100\n')
write('coverage2', 'ID predicted\na 4\nb 1\n')
write('coverage3', 'ID predicted\nb 1\n')
base = ['-E', 'exp', '-d', 'coverage1', 'coverage2', 'coverage3']
assert run(['-n', '1', *base]).stdout == run(base).stdout
assert rows(run(['-n', '2', *base])) == {'a': 3, 'b': 1}
assert rows(run(['-n', '3', *base])) == {'b': 1}
assert rows(run(['-n', '4', *base])) == {}
assert rows(run(['-n', '2', '-w', '1', *base])) == {'a': 3}
assert rows(run(['-n', '2', '-W', '2', *base])) == {'a': 3}
assert set(rows(run(['-n', '2', '-w', '2', '-b', *base]))) == {'a', 'b'}
for count, eligible in [('2', {'a', 'b'}), ('4', set())]:
    run(['-n', count, '-Y', 'coverage_summary', '-T', 'coverage_proto', *base])
    records = Path('coverage_summary').read_text().splitlines()[1:]
    assert {record.split()[0] for record in records} == eligible
    assert set(re.findall(r'id: "([^"\n]+)"', Path('coverage_proto.textproto').read_text())) == eligible
write('invalid_coverage', 'ID predicted\na 2\na nonsense\nb 1\nb 1\n')
assert rows(run(['-E', 'exp', '-d', '-F', '-n', '2', 'invalid_coverage'])) == {'b': 1}
write('repeated_coverage', 'ID predicted\na 2\na 4\n')
assert rows(run(['-E', 'exp', '-d', '-n', '2', 'repeated_coverage'])) == {'a': 3}
write('coverage_class_exp', 'ID class\na A\nb B\nc C\n')
write('coverage_class1', 'ID class\na X\nb B\nc C\n')
write('coverage_class2', 'ID class\nb X\n')
class_base = ['-C', '-E', 'coverage_class_exp', '-d', 'coverage_class1', 'coverage_class2']
assert run(['-n', '1', *class_base]).stdout == run(class_base).stdout
assert run(['-n', '2', *class_base]).stdout.splitlines() == ['ID EXPT N Correct Accuracy', 'b B 2 1 0.5']
assert run(['-n', '3', *class_base]).stdout.splitlines() == ['ID EXPT N Correct Accuracy']
for invalid in ['0', '-1', 'nonsense']:
    run(['-n', invalid, *base], success=False)
# Accumulator precision: identical errors have zero spread, even at large offsets.
for value in ['0.1', '7.1', '1e20']:
    write('constant_predictions', f'ID predicted\na {value}\na {value}\n')
    run(['-E', 'one_exp', '-Y', 'constant_summary', '-T', 'constant_proto', 'constant_predictions'])
    assert float(Path('constant_summary').read_text().splitlines()[1].split()[-1]) == 0.0
    assert 'std_difference: 0' in Path('constant_proto.textproto').read_text()
import struct
first = struct.unpack('f', struct.pack('f', 10000.1))[0]
second = struct.unpack('f', struct.pack('f', 10000.2))[0]
write('offset_predictions', 'ID predicted\na 10000.1\na 10000.2\n')
run(['-E', 'one_exp', '-Y', 'offset_summary', 'offset_predictions'])
stddev = float(Path('offset_summary').read_text().splitlines()[1].split()[-1])
assert abs(stddev - abs(first - second) / (2.0 ** .5)) < 1e-6

# Equal class populations map ties lexically, independently of input order.
write('review_class_pred', 'ID class\na -1\nb 1\n')
expected = ['ID EXPT N Correct Accuracy', 'a A 1 1 1', 'b B 1 1 1']
for data in ['ID class\na A\nb B\n', 'ID class\nb B\na A\n']:
    write('review_class_exp', data)
    actual = run(['-E', 'review_class_exp', '-X', '-1,1', '-d', 'review_class_pred'])
    assert actual.stdout.splitlines() == expected
run(['-E', 'review_class_exp', '-X', '1,1', 'review_class_pred'], success=False)

# Unsupported classification options fail before truncating output files.
write('unsupported_summary', 'keep summary\n')
write('unsupported_proto.textproto', 'keep proto\n')
for mode in [['-C'], ['-c', '0'], ['-X', '-1,1']]:
    for option in [['-w', '1'], ['-W', '1'], ['-Y', 'unsupported_summary'], ['-T', 'unsupported_proto']]:
        run(['-E', 'review_class_exp', *mode, *option, 'review_class_pred'], success=False)
assert Path('unsupported_summary').read_text() == 'keep summary\n'
assert Path('unsupported_proto.textproto').read_text() == 'keep proto\n'
run(['-E', 'one_exp', '-D', 'unsupported_distance', 'one_pred'], success=False)
assert not Path('unsupported_distance').exists()

# Classification -S supplies the same prefix promised by its descriptor header.
write('review_class_pred', 'ID class\na A\nb B\n')
write('review_smiles', 'C a\nCC b\n')
actual = run(['-C', '-E', 'review_class_exp', '-S', 'review_smiles', '-d', 'review_class_pred'])
assert actual.stdout.splitlines() == ['SMILES ID EXPT N Correct Accuracy', 'C a A 1 1 1', 'CC b B 1 1 1']

# Snapshot regex inputs before creating outputs that also match the expression.
write('scan_pred_1', 'ID predicted\na 1\n')
write('scan_pred_2', 'ID predicted\na 3\n')
actual = run(['-E', 'one_exp', '-d', '-Y', 'scan_pred_summary', '-R', '^scan_pred_'])
assert rows(actual) == {'a': 2}
Path('scan_pred_directory').mkdir()
actual = run(['-E', 'one_exp', '-d', '-R', '^scan_pred_(1|2|directory)$'])
assert rows(actual) == {'a': 2}
for name in ['scan_class_1', 'scan_class_2']:
    write(name, Path('review_class_pred').read_text())
actual = run(['-C', '-E', 'review_class_exp', '-d', '-R', '^scan_class_[12]$'])
assert actual.stdout.splitlines() == ['ID EXPT N Correct Accuracy', 'a A 2 2 1', 'b B 2 2 1']

# Fault tolerance applies to records, not missing input lists.
run(['-E', 'one_exp', '-F', 'one_pred', 'F:missing.list'], success=False)
# -k retains the distance and reports it consistently in every continuous output.
write('distance_pred', 'ID predicted distance\na 2 0.2\nb 1 0.5\nc 10 0.9\n')
for fmt in ['def', 'default']:
    result = run(['-E', 'exp', '-d', '-k', fmt, '-Y', 'distance_summary',
                  '-T', 'distance_proto', 'distance_pred'])
    lines = result.stdout.splitlines()
    assert lines[0] == 'ID EXPT AVDIFF NN_DIST', result.stdout
    observed = {f[0]: (float(f[2]), float(f[3])) for f in map(str.split, lines[1:])}
    assert observed == {'a': (2, .2), 'b': (1, .5), 'c': (10, .9)}, observed
    summary = Path('distance_summary').read_text().splitlines()
    assert summary[0].endswith(' NN_DIST'), summary
    assert {f[0]: float(f[-1]) for f in map(str.split, summary[1:])} == {'a': .2, 'b': .5, 'c': .9}
    proto = Path('distance_proto.textproto').read_text()
    assert proto.count('nearest_neighbour_distance:') == 3, proto
result = run(['-E', 'exp', '-k', 'def', 'distance_pred'])
assert ' NN_DIST 0.2' in result.stdout, result.stdout
write('distance_later', 'ID predicted distance\na 4 0.7\n')
result = run(['-E', 'exp', '-d', '-k', 'def', 'distance_pred', 'distance_later'])
assert any(line.split()[0] == 'a' and float(line.split()[-1]) == .7
           for line in result.stdout.splitlines()[1:]), result.stdout
for distance in ['bad', '-0.1', '1.1', 'nan', 'inf']:
    write('invalid_distance', f'ID predicted distance\na 2 {distance}\n')
    run(['-E', 'exp', '-k', 'def', 'invalid_distance'], success=False)
write('missing_distance', 'ID predicted distance\na 2\n')
run(['-E', 'exp', '-k', 'def', 'missing_distance'], success=False)
run(['-E', 'exp', '-k', 'other', 'distance_pred'], success=False)
run(['-E', 'exp', '-k', 'distance_pred'], success=False)
run(['-C', '-E', 'exp', '-k', 'def', 'distance_pred'], success=False)

# /dev/full provides deterministic write failures on systems which have it.
if Path('/dev/full').exists():
    Path('full.textproto').symlink_to('/dev/full')
    for option, destination in [('-Y', '/dev/full'), ('-T', 'full')]:
        run(['-E', 'one_exp', option, destination, 'one_pred'], success=False)
    write('many_exp', 'ID measured\n' + ''.join(f'id{i} 0\n' for i in range(400)))
    write('many_pred', 'ID predicted\n' + ''.join(f'id{i} 1\n' for i in range(400)))
    for option, destination in [('-Y', '/dev/full'), ('-T', 'full')]:
        run(['-E', 'many_exp', option, destination, 'many_pred'], success=False)

print('mispredicted correctness tests passed')
