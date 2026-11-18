#!/usr/bin/env python3
"""Exercise optional InChI grouping with either an enabled or a default build."""
from pathlib import Path
import re
import subprocess
import sys

exe = str(Path(sys.argv[1]).resolve())
source = Path(sys.argv[2])/'in.smi'


def run(*args, success=True):
    result = subprocess.run([exe, *map(str, args)], capture_output=True, text=True)
    if success:
        assert result.returncode == 0, (args, result.returncode, result.stderr)
    else:
        assert result.returncode != 0, (args, result.stdout, result.stderr)
    return result


def groups(text):
    return [line.split(maxsplit=1)[1] for line in text.splitlines()]


def write(name, content):
    path = Path(name)
    path.write_text(content)
    return path


probe = subprocess.run([exe, '-C', '-f', source], capture_output=True, text=True)
if 'does not include InChI' in probe.stderr:
    assert probe.returncode != 0
    assert 'Rebuild with bazel --config=inchi' in probe.stderr
    # Default builds retain normal grouping and fail before creating -C output.
    assert len(groups(run('-f', source).stdout)) == 3
    run('-C', '-S', 'unavailable', source, success=False)
    assert not Path('unavailable.smi').exists()
    print('PASS')
    sys.exit(0)
assert probe.returncode == 0, probe.stderr
assert groups(probe.stdout) == ['pyridone:hydroxypyridine', 'ethanol:ethanol_again']
assert len(groups(run('-f', source).stdout)) == 3
run('-C', '-S', 'normal', source)
assert groups(Path('normal.smi').read_text()) == groups(probe.stdout)
run('-C', '-y', '-S', 'counts', source)
assert groups(Path('counts.smi').read_text()) == ['pyridone 2', 'ethanol 2']
run('-C', '-X', 'num', '-S', 'all_counts', source)
assert groups(Path('all_counts.smi').read_text()) == [
    'pyridone:hydroxypyridine 2', 'ethanol:ethanol_again 2']

bad = run('-C', '-a', '-f', source, success=False)
assert 'incompatible' in bad.stderr

# Previous-file keys must follow the same comparison rule in both modes.
previous = write('previous.smi', 'OC1=NC=CC=C1 previous\n')
for mode in [[], ['-f']]:
    run('-C', *mode, '-p', previous, '-p', 'suffix=_old', '-S', 'previous_output', source)
    assert groups(Path('previous_output.smi').read_text()) == [
        'previous_old:pyridone:hydroxypyridine', 'ethanol:ethanol_again']

# Stereo suppression changes the comparison only, retaining the representative.
stereo = write('stereo.smi', 'F[C@H](Cl)Br r\nF[C@@H](Cl)Br s\n'
               'F/C=C/F trans\nF/C=C\\F cis\n')
for mode in [[], ['-f']]:
    for flags, expected in [([], ['r', 's', 'trans', 'cis']),
                            (['-c'], ['r:s', 'trans', 'cis']),
                            (['-x'], ['r', 's', 'trans:cis']),
                            (['-c', '-x'], ['r:s', 'trans:cis']),
                            (['-c', '-Z'], ['r:s', 'trans:cis'])]:
        run('-C', *mode, *flags, '-S', 'stereo_output', stereo)
        output = Path('stereo_output.smi').read_text()
        assert groups(output) == expected, (mode, flags, output)
        assert '@' in output and '/' in output, output

# Existing isotope and largest-fragment preprocessing applies before InChI.
parents = write('parents.smi', '[13CH3]CO labelled\nCCO parent\nCCO.[Na+] salt\n')
for flags, expected in [([], ['labelled', 'parent', 'salt']),
                        (['-I'], ['labelled:parent', 'salt']),
                        (['-I', '-l'], ['labelled:parent:salt'])]:
    assert groups(run('-C', '-f', *flags, parents).stdout) == expected

# Proto keys must be InChIKeys too; molecule/name reporting remains unchanged.
run('-C', '-X', 'proto=groups.textproto', source)
proto = Path('groups.textproto').read_text()
assert len(re.findall(r'key: "[A-Z]{14}-[A-Z]{10}-[A-Z]"', proto)) == 2, proto
assert 'id: "pyridone"' in proto and 'id: "hydroxypyridine"' in proto

# Unsupported structures must fail, not accumulate under an empty key.
invalid = write('unsupported.smi', '* wildcard\n')
for mode in [[], ['-f']]:
    error = run('-C', *mode, '-S', 'invalid_output', invalid, success=False)
    assert "Cannot generate InChI comparison key for 'wildcard'" in error.stderr, error.stderr
error = run('-C', '-f', '-p', invalid, source, success=False)
assert "Cannot generate InChI comparison key for 'wildcard'" in error.stderr
print('PASS')
