#!/usr/bin/env python3
"""Regressions for counts, seed distances, weights, squeezing and stable ties."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())
fixture = Path(__file__).resolve().parents[2] / 'gfp_leader_standard/reassignment/in/fingerprints.gfp'
records = [r + '|\n' for r in fixture.read_text().split('|\n') if r.strip()]

def parse(text):
    result = []
    for block in text.split('|\n'):
        names = re.findall(r'^PCN<([^>]+)>', block, re.M)
        if names:
            result.append((names[0], re.search(r'^DIST<([^>]+)>', block, re.M)[1]))
    return result

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    def write(name, text):
        path = tmp / name
        path.write_text(text)
        return path
    def renamed(record, name):
        return re.sub(r'PCN<[^>]+>', f'PCN<{name}>', record)
    def run(path, *args, success=True):
        r = subprocess.run([exe, *([] if '-h' in args else ['-h', '1']), *args, str(path)], capture_output=True,
                           text=True, timeout=30)
        assert (r.returncode == 0) == success, r.stderr
        if not success:
            assert not r.stdout
        return r

    same = write('same.gfp', ''.join(renamed(records[0], f'm{i}') for i in range(5)))
    ones = write('ones.txt', 'ID Weight\n' + ''.join(f'm{i} 1\n' for i in range(5)))
    for weighted in ([], ['-p', f'FILE={ones}']):
        for n in (1, 2, 4, 5):
            selected = parse(run(same, '-n', str(n), *weighted).stdout)
            assert [name for name, _ in selected] == [f'm{i}' for i in range(n)]
            assert [float(d) for _, d in selected] == [1] + [0] * (n - 1)
        assert len(parse(run(same, '-n', '5', '-t', '.01', *weighted).stdout)) == 1

    pool = write('pool.gfp', ''.join(renamed(r, f'm{i}') for i, r in enumerate(records[:5])))
    uniform = run(pool, '-n', '5').stdout
    weighted = run(pool, '-n', '5', '-p', f'FILE={ones}').stdout
    assert parse(uniform) == parse(weighted)
    assert len(parse(run(pool, '-t', '.99', '-p', f'FILE={ones}').stdout)) == 1

    # All -N input forms supply the same initial distances, including names with
    # extra columns when looking up a plain id from a numeric data file.
    seeds = [0, .1, .2, .3, .4]
    column = write('column.gfp', ''.join(renamed(r, f'm{i} {seeds[i]}')
                                        for i, r in enumerate(records[:5])))
    tagged = write('tagged.gfp', column.read_text().replace('|\n', 'SEED<placeholder>\n|\n'))
    text = tagged.read_text()
    for d in seeds:
        text = text.replace('SEED<placeholder>', f'SEED<{d}>', 1)
    tagged.write_text(text)
    seedfile = write('seed.txt', 'ID Distance\n' + ''.join(f'm{i} {d}\n' for i, d in enumerate(seeds)))
    first = []
    for args in (['-N', 'COL=2'], ['-N', 'SEED'], ['-N', f'FILE={seedfile}']):
        text = run(tagged, '-n', '1', *args).stdout
        first.append(parse(text))
    assert first[0] == first[1] == first[2] == [('m4 0.4', '0.4')]
    prior_last = write('prior_last.gfp', records[4])
    combined = run(tagged, '-n', '1', '-N', 'SEED', '-A', str(prior_last)).stdout
    assert parse(combined) == [('m3 0.3', '0.3')]
    seed_and_weight = write('seed_weight.gfp', ''.join(
        renamed(r, f'm{i} {seeds[i]} {i+1}') for i, r in enumerate(records[:5])))
    weighted_seed = run(seed_and_weight, '-n', '1', '-N', 'COL=2', '-p', 'COL=3').stdout
    assert parse(weighted_seed) == [('m4 0.4 5', '2')]
    zero = run(same, '-n', '1', '-A', str(same), '-p', f'FILE={ones}').stdout
    assert parse(zero)[0][1] == '0'

    # COL=2 and COL=3 must agree with equivalent full-file columns.
    scales = [1, .5, 2, .25, 4]
    columns = write('columns.gfp', ''.join(renamed(r, f'm{i} {scales[i]} {scales[4-i]}')
                                         for i, r in enumerate(records[:5])))
    scale_file = write('scales.txt', 'ID A B\n' + ''.join(f'm{i} {w} {scales[4-i]}\n'
                                                       for i, w in enumerate(scales)))
    for col in (2, 3):
        actual = run(columns, '-p', f'COL={col}').stdout
        expected = run(columns, '-p', f'FILE={scale_file}', '-p', f'FCOL={col}').stdout
        assert actual == expected
    assert run(columns, '-p', 'COL=2').stdout != run(columns, '-p', 'COL=3').stdout
    tabs = write('scales.tsv', scale_file.read_text().replace(' ', '\t'))
    assert run(columns, '-p', f'FILE={tabs}').stdout == run(
        columns, '-p', f'FILE={scale_file}').stdout

    # Removing a component must not turn identical structures into distant ones.
    def empty_mk2(record):
        def replace(match):
            fields = match[1].split(';')
            fields[0] = '.' * (len(fields[0]) - 1) + fields[0][-1]
            fields[2] = fields[4] = '0'
            return 'FPMK2<' + ';'.join(fields) + '>'
        return re.sub(r'FPMK2<([^>]+)>', replace, record)
    empty = write('empty.gfp', empty_mk2(same.read_text()))
    assert [float(d) for _, d in parse(run(empty, '-n', '2').stdout)] == [1, 0]

    # Squeezing must preserve original ids, neighbour ids, scales and distances.
    large = write('large.gfp', ''.join(renamed(records[i % 12], f'm{i} {scales[i % 5]}')
                                      for i in range(1104)))
    prior = write('prior.gfp', records[0])
    for weighted in ([], ['-p', 'COL=2']):
        for previous in ([], ['-A', str(prior)]):
            baseline = run(large, '-n', '30', '-q', '0', *weighted, *previous).stdout
            for threads in ('1', '4'):
                squeezed = run(large, '-n', '30', '-q', '1', '-h', threads,
                               *weighted, *previous).stdout
                assert squeezed == baseline
            assert len(parse(baseline)) == 30

    malformed = write('bad.gfp', re.sub(r'FPIW<[^>]+>', 'FPIW<bad>', records[0]) + records[1])
    run(malformed, success=False)
    run(tmp / 'missing.gfp', success=False)
    run(write('empty.gfp', ''), success=False)
    run(write('not_tdt.gfp', 'not a fingerprint file\n'), success=False)
    for text in ('m0 1\n' * 5, 'ID W\n' + ''.join(f'm{i} nan\n' for i in range(5)),
                 'ID W\n' + ''.join(f'm{i} inf\n' for i in range(5))):
        weights = write('bad_weights.txt', text)
        run(pool, '-p', f'FILE={weights}', success=False)
    for token in ('nan', 'inf', '-1', '2'):
        invalid = write('bad_seed.gfp', ''.join(renamed(r, f'm{i} {token}')
                                              for i, r in enumerate(records[:5])))
        run(invalid, '-N', 'COL=2', success=False)
    run(pool, '-h', '0', success=False)
print('gfp_spread_standard correctness checks passed')
