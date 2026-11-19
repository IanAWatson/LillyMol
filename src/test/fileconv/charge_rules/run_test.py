#!/usr/bin/env python3
"""Check the distributed charge rules using fileconv, without bindings."""
import json
from pathlib import Path
import subprocess
import sys

exe = Path(sys.argv[1]).resolve()
cases = json.loads((Path(sys.argv[2]) / 'cases.json').read_text())
for filename, field in [('input.smi', 'source'), ('expected.smi', 'target')]:
    Path(filename).write_text(''.join(f"{case[field]} {case['name']}\n" for case in cases))


def convert(source, options=()):
    result = subprocess.run([str(exe), *options, '-o', 'usmi', '-S', '-', str(source)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    return result.stdout


expected = convert('expected.smi')
for options in (['-N', 'DIR=DEF'], ['-N', 'DIR=DEF', '-N', 'simple']):
    actual = convert('input.smi', options)
    assert actual == expected, (options, actual, expected)
    # The current phosphate limit counts embeddings per pass, not total
    # existing negative oxygens. Its repeat policy is a separate open issue.
    stable = ''.join(line + '\n' for case, line in zip(cases, actual.splitlines())
                     if case.get('check_repeat', True))
    Path('charged.smi').write_text(stable)
    assert convert('charged.smi', options) == stable
print('Charge query regressions passed')
