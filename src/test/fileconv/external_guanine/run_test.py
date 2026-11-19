#!/usr/bin/env python3
"""Exercise the distributed external rules using CLI tools, without bindings."""
import json
import os
from pathlib import Path
import subprocess
import sys

exe = Path(sys.argv[1]).resolve()
indir = Path(sys.argv[2])
root = Path(os.environ.get('LILLYMOL_HOME', Path(__file__).resolve().parents[4]))
rules = root / 'data/standardisation/guanine.textproto'
cases = json.loads((indir / 'cases.json').read_text())
Path('input.smi').write_text(''.join(f'{source} {name}\n' for name, source, _ in cases))
Path('expected.smi').write_text(''.join(f'{target} {name}\n' for name, _, target in cases))


def convert(source, options=()):
    result = subprocess.run([str(exe), *options, '-o', 'usmi', '-S', '-', str(source)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return result.stdout


external = ['-g', f'EXT:{rules}']
# Actual and target molecules are independently canonicalised. The all case
# also checks that unsupported external cases retain built-in behaviour.
for prefix in ([], ['-g', 'all']):
    # Built-in imidazole processing can choose the other N-H position.
    # Supported guanines must finish in the explicitly requested external
    # target, while unsupported cases retain the built-in result.
    canonical_targets = convert('expected.smi').splitlines()
    built_in_targets = convert('expected.smi', prefix).splitlines()
    supported = {'keto', 'enol', '7H', 'imino', 'N9_methyl', 'N7_methyl',
                 'isotope', 'remote_chiral', 'multiple'}
    expected = ''.join((canonical if name in supported else built_in) + '\n'
                       for (name, _, _), canonical, built_in in
                       zip(cases, canonical_targets, built_in_targets))
    actual = convert('input.smi', [*prefix, *external])
    assert actual == expected, (prefix, actual, expected)
    Path('roundtrip.smi').write_text(actual)
    assert convert('roundtrip.smi', [*prefix, *external]) == actual

# The external rules alone exclude anionic oxygen. With -g all the built-in
# rules may neutralise it first, making it eligible, so check this separately.
Path('charged.smi').write_text('NC1=NC2=C(N=CN2)C([O-])=N1 oxygen_anion\n')
assert convert('charged.smi', external) == convert('charged.smi')

# Neither the CLI return status nor output equality alone detects every
# tstandardise failure. Inspect its diagnostics and failure file too.
checker = exe.with_name('tstandardise')
for name, prefix in [('external', []), ('all', ['-g', 'all'])]:
    failures = Path(f'{name}.fail.smi')
    result = subprocess.run([str(checker), *prefix, *external, '-n', '100',
                             '-s', '7341', '-v', '-F', str(failures), 'input.smi'],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert failures.read_text() == '', result.stderr
    assert '0 molecules with failures\n' in result.stderr, result.stderr
    assert '0 valence errors introduced\n' in result.stderr, result.stderr
print('External guanine tests passed')
