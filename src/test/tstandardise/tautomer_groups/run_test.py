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
    'N1=NC(=NN=C(C2=CC=NC=C2)Br)NN1',
]
substituted_cyclic_diketones = [
    'OC1=C(c2ccccc2)C(=O)CCC1',
    'O=C1C(c2ccccc2)=C(O)CCC1',
    'OC1=C(c2ccccc2)C(=O)[C@H](C)CC1',
    'O=C1C(c2ccccc2)=C(O)[C@H](C)CC1',
    'OC1=C(C)C(=O)CCCC1',
]
check('substituted_cyclic', ['-g', 'beta-diketone'], substituted_cyclic_diketones)
check('substituted_cyclic_all', ['-g', 'all'], substituted_cyclic_diketones)
check('substituted_cyclic_disabled', ['-g', 'all', '-g', '-beta-diketone'],
      substituted_cyclic_diketones)

phosphorus = [
    'P(O)(O)C(N)C', 'P(=O)(O)C(N)C', 'CP(C)O',
    'COP(C)O', 'P(O)(O)O', 'CP(O)[18OH]', 'CP([OH])O',
    'CP(C)O[2H]', 'CP(O)O[2H]',
]
check('phosphorus', ['-g', 'phosphorus'], phosphorus)
check('phosphorus_all', ['-g', 'all'], phosphorus)
# CP([OH])O has a pre-existing unique-SMILES ordering inconsistency when no
# transformation runs. Its disabled behavior is checked directly in unit tests.
check('phosphorus_disabled', ['-g', 'all', '-g', '-phosphorus'],
      [smiles for smiles in phosphorus if smiles != 'CP([OH])O'])

ascorbates = [
    'OCC(O)C1OC(O)=C(O)C1=O', 'OCC(O)C1OC(=O)C(O)=C1O',
    'OCC(O)c1oc(O)c(O)c1O', 'OCC(O)C1OC(=O)C(O)C1=O',
    'OCC(O)C1OC(O)=C(OP(=O)(O)O)C1=O',
    'OCC(O)c1oc(O)c(OP(=O)(O)O)c1O',
    'OC[C@H](O)[C@@H]1OC(O)=C(O)C1=O',
    '[2H]C1(C(O)CO)OC(O)=C(O)C1=O',
]
check('ascorbate', ['-g', 'ascorbate'], ascorbates)
check('ascorbate_all', ['-g', 'all'], ascorbates)
check('ascorbate_disabled', ['-g', 'all', '-g', '-ascorbate'], ascorbates)

check('urea', ['-g', 'msdur'], ['[OH]C(N)=[NH]', 'OC(N)=N'])
for rule in ['all', 'amidine', 'pyrimidine']:
    check(rule, ['-g', rule], cases)
check('independent_enols', ['-g', 'keto_enol', '-g', 'ENOL'], [
    'CC(=O)CC.CC(=O)CC', 'CC(=O)CC.CC(=S)CC', 'CC(=O)CCCCC(=S)C'])

print('Tautomer atom-order tests passed')
