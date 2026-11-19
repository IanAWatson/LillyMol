#!/usr/bin/env python3
"""Test the distributed hypoxanthine templates without Python bindings."""
import json
import os
import platform
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile

exe = Path(sys.argv[1]).resolve()
indir = Path(sys.argv[2])
root = Path(os.environ.get('LILLYMOL_HOME', Path(__file__).resolve().parents[4]))
rules = root / 'data/standardisation/hypoxanthine.textproto'
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

# Charges are excluded by EXT itself, but built-in neutralisation may make
# them eligible when -g all is used. Check the exclusion independently.
Path('charged.smi').write_text('[O-]C1=NC=NC2=C1N=CN2 charged\n')
assert convert('charged.smi', external) == convert('charged.smi')

# Guanine and hypoxanthine specifications can be loaded together in either
# order. Their distinct C2 substitution requirements keep their scopes apart.
guanine = ['-g', f'EXT:{root / "data/standardisation/guanine.textproto"}']
Path('combined.smi').write_text('OC1=NC=NC2=C1N=CN2 hypoxanthine\n'
                              'NC1=NC2=C(N=CN2)C(O)=N1 guanine\n')
a = convert('combined.smi', ['-g', 'all', *guanine, *external])
b = convert('combined.smi', ['-g', 'all', *external, *guanine])
assert a == b, (a, b)
# Exercise the convenience wrapper using the tested executable in a temporary
# installation. Spaces in both the installation and input paths must survive.
with tempfile.TemporaryDirectory(prefix='lillymol standardise ') as directory:
    home = Path(directory)
    bindir = home / 'bin' / platform.system()
    bindir.mkdir(parents=True)
    (bindir / 'fileconv').symlink_to(exe)
    (home / 'data').symlink_to(root / 'data', target_is_directory=True)
    scripts = home / 'contrib/bin'
    scripts.mkdir(parents=True)
    wrapper = scripts / 'standardise_all.sh'
    shutil.copy2(root / 'contrib/bin/standardise_all.sh', wrapper)
    spaced_input = home / 'input molecules.smi'
    spaced_input.write_text(Path('combined.smi').read_text())
    for infer_home in (False, True):
        env = os.environ.copy()
        if infer_home:
            env.pop('LILLYMOL_HOME', None)
        else:
            env['LILLYMOL_HOME'] = str(home)
        result = subprocess.run([str(wrapper), '-o', 'usmi', '-S', '-', str(spaced_input)],
                                env=env, capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stderr
        assert result.stdout == a, (infer_home, result.stdout, a)
    result = subprocess.run([str(wrapper), '--help'], capture_output=True, text=True)
    assert result.returncode == 0 and 'EXT' in result.stdout, result.stderr
    result = subprocess.run([str(wrapper), str(home / 'missing.smi')], env=env,
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stderr

print('External hypoxanthine tests passed')
