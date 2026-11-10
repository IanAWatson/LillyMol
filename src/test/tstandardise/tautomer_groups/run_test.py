#!/usr/bin/env python3
"""Check random atom order and repeated processing of interacting tautomer rules."""
from pathlib import Path
import subprocess
import sys


def check(name, options, smiles):
    source = Path(name + '.smi')
    source.write_text(''.join(f'{smi} molecule_{i}\n' for i, smi in enumerate(smiles)))
    failures = Path(name + '.fail.smi')
    result = subprocess.run(
        [sys.argv[1], *options, '-n', '100', '-v', '-F', str(failures), str(source)],
        capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    # tstandardise can return success even when it detects inconsistencies.
    assert failures.read_text() == '', result.stderr
    assert '0 molecules with failures\n' in result.stderr, result.stderr
    assert '0 valence errors introduced\n' in result.stderr, result.stderr


cases = [
    'CC(N)=NC', 'CC(=N)NC', 'C(N)=NC', 'CC(N)=NC1CCCCC1',
    'CC(N)=Nc1ccccc1', 'CC(N)=NC(=N)N', 'C/C(N)=N/C', 'CC(N)=[NH+]C',
    'N1=C(O)C=C(O)N=C1O', 'N1C(=O)CC(=O)NC1=O',
    'N1C(=O)C=C(O)N=C1O', 'N1=C(S)C=C(O)N=C1O',
    'N1=C(S)C=C(S)N=C1S', 'CN1C(O)=CC(=O)N(C)C1=O',
    'N1=C(O)C(C)(C)C(O)=NC1=O', 'N1C(O)CC(O)NC1O',
    'N1=C(O)C=C(O)N=C1O.O=c1cc[nH]c(C)n1',
]
for rule in ['all', 'amidine', 'pyrimidine']:
    check(rule, ['-g', rule], cases)
check('independent_enols', ['-g', 'keto_enol', '-g', 'ENOL'], [
    'CC(=O)CC.CC(=O)CC', 'CC(=O)CC.CC(=S)CC', 'CC(=O)CCCCC(=S)C'])

print('Tautomer atom-order tests passed')
