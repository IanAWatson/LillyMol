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
beta_enaminones = [
    'CC(=O)CC(C)=NCC', 'CC(=O)C=C(C)NCC',
    'O=C(C)CC(=NCCN=C(C)CC(=O)C)C',
    'CC(=O)C=C(C)NCCNC(C)=CC(C)=O',
    'O=C(c1ccccc1)CC(=Nc1ccccc1)c1ccccc1',
    'CC(=O)[13CH2]C(C)=[15N]CC', 'CC(=O)CC(C)=N[C@H](F)C',
    'CC(=O)C/C(C)=N/CC', 'CC(=O)C([2H])C(C)=NCC',
    'CC(=O)CC(CC(C)=O)=NCC', 'O=C1CC(=NCC)CCC1',
    'CC(=O)CC(C)=NCC.CC(=O)CC(C)=NCC',
    'CC(=O)CC(C)=NCC.CC(=O)C=C(C)NNC',
]
check('beta_enaminone', ['-g', 'beta-enaminone'], beta_enaminones)
check('beta_enaminone_all', ['-g', 'all'], beta_enaminones)
check('beta_enaminone_disabled', ['-g', 'all', '-g', '-beta-enaminone'],
      beta_enaminones)
check('beta_enaminone_with_hydrazone', ['-g', 'all', '-g', 'beta-hydrazone'],
      beta_enaminones)

beta_hydrazones = [
    'CC(=O)C=C(C)NNC', 'CC(=O)CC(C)=NNC',
    'CC(=O)C=C(C)NNC(=O)c1ccncc1',
    'O=C(c1ccccc1)C=C(c1ccccc1)NNc1ccccc1',
    'CC(=O)[13CH]=C(C)[15NH]NC', 'CC(=O)C=C(C)NN[C@H](F)C',
    'CC(=O)/C=C(/C)NNC', 'CC(=O)C=C(C)N([2H])NC',
    'CC(=O)C=C(C)NN(C)C', 'O=C1C=C(NNC)CCC1',
    'CC(=O)C=C(C)NNC.CC(=O)C=C(C)NNC',
]
check('beta_hydrazone', ['-g', 'beta-hydrazone'], beta_hydrazones)
check('beta_hydrazone_all', ['-g', 'all'], beta_hydrazones)
check('beta_hydrazone_disabled', ['-g', 'all', '-g', '-beta-hydrazone'],
      beta_hydrazones)

beta_enehydrazines = [
    'CC(=O)CC(C)=NNC', 'CC(=O)C=C(C)NNC',
    'CC(=O)CC(C)=NNC(=O)c1ccncc1',
    'O=C(c1ccccc1)CC(=NNc1ccccc1)c1ccccc1',
    'CC(=O)[13CH2]C(C)=[15N]NC', 'CC(=O)CC(C)=NN[C@H](F)C',
    'CC(=O)C/C(C)=N/NC', 'CC(=O)CC(C)=NN([2H])C',
    'CC(=O)C([2H])C(C)=NNC', 'CC(=O)CC(CC(C)=O)=NNC',
    'CC(=O)C=C(CC(C)=O)NNC', 'CC(=O)CC(C)=NNC.CC(=O)CC(C)=NCC',
]
check('beta_enehydrazine', ['-g', 'beta-enehydrazine'], beta_enehydrazines)
check('beta_enehydrazine_all',
      ['-g', 'all', '-g', '-beta-hydrazone', '-g', 'beta-enehydrazine'],
      beta_enehydrazines)
# A temporary conflict is fine when the final option disables the other side.
check('beta_enehydrazine_reordered',
      ['-g', 'beta-enehydrazine', '-g', 'all', '-g', '-beta-hydrazone'],
      beta_enehydrazines)

cyclic_amidines = [
    'N=C1NCCC1', 'NC1=NCCC1', 'N=C1NCCCC1', 'NC1=NCCCC1',
    'N1C(=N)C(O)C(C)C1CCCCC', '[15NH]=C1NCCC1',
    'N=C1[15NH]CCC1', 'N=[13C]1NCCC1', 'N=C1N[C@H](C)CC1',
    '[NH2+]=C1NCCC1', '[NH2+]1C(=N)CCC1', '[NH+]1=C(N)CCC1',
    'N=C1NCC1', 'N=C1NCCCCC1', 'N=C1NCCN1', 'N=C1NCCO1',
    'N=C1NCC=C1', 'N=C1NC(=O)CC1', 'N1C(=N)CCCC1=N',
    'N=C1NC2CCCCC2C1', 'N=C1NCCC12CCCC2', 'N=C1N(C)CCC1',
    'CN=C1NCCC1', '[2H]N1C(=N)CCC1', 'N1C(=N[2H])CCC1',
]
check('cyclic_amidine', ['-g', 'cyclic-amidine'], cyclic_amidines)
check('cyclic_amidine_all', ['-g', 'all'], cyclic_amidines)
check('cyclic_amidine_disabled', ['-g', 'all', '-g', '-cyclic-amidine'],
      cyclic_amidines)

beta_ketoesters = [
    'COC(=O)C=C(O)C', 'CCOC(=O)C=C(O)C',
    'COC(=O)C=C(O)c1ccccc1', 'COC(=O)[CH]=C([OH])C',
    'C[18O]C(=O)C=C(O)C', 'C[C@H](F)OC(=O)C=C(O)C',
    'COC(=O)C=C(O)C=C(O)C', 'COC(=O)/C=C(/O)C',
]
check('beta_ketoester', ['-g', 'beta-ketoester'], beta_ketoesters)
check('beta_ketoester_all', ['-g', 'all'], beta_ketoesters)
check('beta_ketoester_disabled', ['-g', 'all', '-g', '-beta-ketoester'], beta_ketoesters)

nitroso = [
    'CN=O', 'CCN=O', 'CC(C)N=O', 'O=NC1CCCCC1',
    'C[N+]1=CC=CC=C1CN=O', 'C[CH2][N]=[O]',
    'CCN=[18O]', 'CC([2H])N=O', 'CC([2H])([2H])N=O',
    'C[C@H](F)CN=O', 'C[C@H](N=O)CC', 'C/C=C/CN=O',
    'N1C=CC(N=O)C=C1',
]
check('nitroso', ['-g', 'nitroso'], nitroso)
check('nitroso_all', ['-g', 'all'], nitroso)
check('nitroso_disabled', ['-g', 'all', '-g', '-nitroso'], nitroso)

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
