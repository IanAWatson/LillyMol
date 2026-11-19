#!/usr/bin/env python3
"""Test the distributed allopurinol templates without Python bindings."""
import json
import os
import platform
from pathlib import Path
import subprocess
import sys
import tempfile

exe = Path(sys.argv[1]).resolve()
indir = Path(sys.argv[2])
root = Path(os.environ.get('LILLYMOL_HOME', Path(__file__).resolve().parents[4]))
rules = root / 'data/standardisation/allopurinol.textproto'
cases = json.loads((indir / 'cases.json').read_text())
for filename, field in [('input.smi', 'source'), ('expected.smi', 'target')]:
    Path(filename).write_text(''.join(f"{case[field]} {case['name']}\n" for case in cases))


def convert(source, options=()):
    result = subprocess.run([str(exe), *options, '-o', 'usmi', '-S', '-', str(source)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return result.stdout


external = ['-g', f'EXT:{rules}']
canonical_targets = convert('expected.smi').splitlines()
for label, prefix in [('external', []), ('all', ['-g', 'all']),
                      ('enol', ['-g', 'all', '-g', 'ENOL'])]:
    built_in_targets = convert('expected.smi', prefix).splitlines()
    expected = ''.join((canonical if case['supported'] else built_in) + '\n'
                       for case, canonical, built_in in
                       zip(cases, canonical_targets, built_in_targets))
    actual = convert('input.smi', [*prefix, *external])
    assert actual == expected, (label, actual, expected)
    Path('roundtrip.smi').write_text(actual)
    assert convert('roundtrip.smi', [*prefix, *external]) == actual
    failures = Path(f'{label}.fail.smi')
    result = subprocess.run([str(exe.with_name('tstandardise')), *prefix, *external,
                             '-n', '100', '-s', '7341', '-v', '-F', str(failures),
                             'input.smi'], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert failures.read_text() == '', result.stderr
    assert '0 molecules with failures\n' in result.stderr, result.stderr
    assert '0 valence errors introduced\n' in result.stderr, result.stderr

# The guards concern the input to EXT, not the original record: built-in
# processing can first make an unsupported form eligible. Check them alone.
guards = json.loads((indir / 'guarded.json').read_text())
Path('guarded.smi').write_text(''.join(f'{smi} {name}\n' for name, smi in guards))
assert convert('guarded.smi', external) == convert('guarded.smi')
Path('charged.smi').write_text('[O-]C1=NC=NC2NN=CC1=2 charged\n')
assert convert('charged.smi', external) == convert('charged.smi')

# Ensure the convenience wrapper really enables this new specification,
# not just the earlier guanine and hypoxanthine rules.
with tempfile.TemporaryDirectory(prefix='lillymol allopurinol ') as directory:
    home = Path(directory)
    bindir = home / 'bin' / platform.system()
    bindir.mkdir(parents=True)
    (bindir / 'fileconv').symlink_to(exe)
    (home / 'data').symlink_to(root / 'data', target_is_directory=True)
    source = home / 'molecules.smi'
    source.write_text('O=C1NC=NC2=NNC=C12 allopurinol\n'
                      'O=C1NC(=O)NC2=NNC=C12 oxypurinol\n')
    env = os.environ.copy()
    env['LILLYMOL_HOME'] = str(home)
    wrapper = root / 'contrib/bin/standardise_all.sh'
    result = subprocess.run([str(wrapper), '-o', 'usmi', '-S', '-', str(source)],
                            env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert result.stdout == convert(source, ['-g', 'all', *external])
print('External allopurinol tests passed')
