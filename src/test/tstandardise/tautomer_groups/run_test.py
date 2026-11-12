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
phosphorus_oxo = [
    'C[P+](C)(C)O', 'CO[P+](N)(N)O', 'C[P+](O)([18OH])OC',
    'C[P+](O)([OH])OC', 'C[P+](O)(O[2H])OC',
    'C[P@+](O)(N)OC', 'C[P+](C)(C)O[2H]', 'C[P+](C)(C)C',
    'C[P+](C)(C)O.C[P+](C)(C)O.C[NH2+]C',
    '[P+]1(O)(NCC(Cl)CO1)N(CCCl)CCCl',
    'C[P+](O)(N)[C@H](F)C',
]
check('phosphorus_oxo', ['-g', 'phosphorus-oxo'], phosphorus_oxo)
check('phosphorus_oxo_all', ['-g', 'all'], phosphorus_oxo)
check('phosphorus_oxo_enol', ['-g', 'all', '-g', 'ENOL'], phosphorus_oxo)
check('phosphorus_oxo_disabled', ['-g', 'all', '-g', '-phosphorus-oxo'],
      phosphorus_oxo)

oxo_ring_guanidines = [
    'O=C1NC(=N)N(C)C1',
    'O=C1NC(=N)NC1',
    'O=C1NC(N)=NC1',
    'O=C1NC(=N)N(c2ccccc2)C1',
    'O=C1NC(=N)N(C)C(C)1',
    'O=C1NC(=N)N(C)[C@H]1C',
    'O=C1[15NH]C(=[15NH])N(C)[13CH2]1',
    'O=C1N=C(N)N(C)C1',
    'O=C1N=C(N)NC1',
    'O=C1NC(=N)NCC1',
    'S=C1NC(=N)NC1',
    'O=C1NC(=N)NNC1',
    'O=C1N=C(N)N(C(=O)C)C1',
    'O=C1N=C(N)N(O)C1',
    'O=C1NC(=N)NC1=O',
    'O=C1NC(=N[2H])NC1',
    'O=C1N([2H])C(=N)NC1',
]
check('oxo_ring_guanidine', ['-g', 'Rguan'], oxo_ring_guanidines)
check('oxo_ring_guanidine_all', ['-g', 'all'], oxo_ring_guanidines)
check('oxo_ring_guanidine_disabled', ['-g', 'all', '-g', '-Rguan'],
      oxo_ring_guanidines)

local_ring_guanidines = [
    'N1C(=Nc2ccccc2)NCC1',
    'N1C(=NC)NCC1',
    'N1C(=Nc2ccccc2)NCCC1',
    'N1C(=Nc2ccccc2)NC(C)C1',
    'N1C(Nc2ccccc2)=NC(C)C1',
    'N1=C(Nc2ccccc2)NC(C)C1',
    'N1C(=Nc2ccccc2)NC(C)CC1',
    'N1C(=Nc2ccccc2)[15NH]CC1',
    'N1C(=Nc2ccccc2)NC[13CH2]1',
    '[NH]1C(=Nc2ccccc2)NCC1',
    'N1C(=N[C@H](F)C)NCC1',
    'N1C(=Nc2ccccc2)NC(C)C1C',
    'N1C(=Nc2ccccc2)NC(C)CC1CC',
    'N1C(=Nc2ccccc2)N[C@H](C)C1',
    'N1C(=NC(=O)C)NCC1',
    'N1C(=N)NCC1',
    'N1C(=Nc2ccccc2)NCCCC1',
    'O=C1NC(=Nc2ccccc2)NC1',
    'N1C(=Nc2ccccc2)NNC1',
]
check('local_ring_guanidine', ['-g', 'Rguan'], local_ring_guanidines)
check('local_ring_guanidine_all', ['-g', 'all'], local_ring_guanidines)
check('local_ring_guanidine_disabled', ['-g', 'all', '-g', '-Rguan'],
      local_ring_guanidines)

alpha_ketoamides = [
    'NC(=O)C(O)=C',
    'CNC(=O)C(O)=CC',
    'CN(C)C(=O)C(O)=Cc1ccccc1',
    'O=C(N1CCOCC1)C(O)=CC',
    'c1ccccc1NC(=O)C(O)=CC',
    'NC(=O)C([18OH])=[13CH]C',
    'C[C@H](F)NC(=O)C(O)=CC',
    'NC(=O)C(=O)CC',
    'ONC(=O)C(O)=CC',
    'NNC(=O)C(O)=CC',
    'CC(=O)NC(=O)C(O)=CC',
    'CC(=S)NC(=O)C(O)=CC',
    'CNC(=N)C(O)=CC',
    'C[NH2+]C(=O)C(O)=CC',
    'COC(=O)C(O)=CC',
    'OC(=O)C(O)=CC',
    'CSC(=O)C(O)=CC',
    'NC(=O)C(O)=CC=C',
    'NC(=O)C(O)=CC(=O)N',
    'NC(=O)C(O)=CC(=O)C',
    'NC(=O)C(O)=C(C)C',
    'NC(=O)C1=C(O)CCC1',
    'O=C1NC(O)=CCC1',
    'NC(=O)/C(O)=C/C',
    'NC(=O)C(O[2H])=CC',
    'NC(=O)C(O)=C([2H])C',
    'NC(=O)C([O-])=CC',
]
check('alpha_ketoamide', ['-g', 'alpha-ketoamide'], alpha_ketoamides)
check('alpha_ketoamide_all', ['-g', 'all'], alpha_ketoamides)
check('alpha_ketoamide_disabled', ['-g', 'all', '-g', '-alpha-ketoamide'],
      alpha_ketoamides)
check('alpha_ketoamide_enol_preference', ['-g', 'all', '-g', 'ENOL'],
      alpha_ketoamides)

alpha_ketoesters = [
    'COC(=O)C(O)=C',
    'CCOC(=O)C(O)=CC',
    'COC(=O)C(O)=Cc1ccccc1',
    'c1ccccc1OC(=O)C(O)=CC',
    'COC(=O)C([18OH])=[13CH]C',
    'C[C@H](F)COC(=O)C(O)=CC',
    'COC(=O)C(=O)CC',
    'OC(=O)C(O)=CC',
    'NC(=O)C(O)=CC',
    'CC(=O)OC(=O)C(O)=CC',
    'C=COC(=O)C(O)=CC',
    'COC(=O)C(O)=CC=C',
    'COC(=O)C(O)=CC(=O)N',
    'COC(=O)C(O)=CC(=O)C',
    'COC(=O)C(O)=C(C)C',
    'COC(=O)C1=C(O)CCC1',
    'O=C1OC(O)=CCC1',
    'COC(=O)/C(O)=C/C',
    'COC(=O)C(O[2H])=CC',
    'COC(=O)C(O)=C([2H])C',
    'COC(=O)C([O-])=CC',
]
check('alpha_ketoester', ['-g', 'alpha-ketoester'], alpha_ketoesters)
check('alpha_ketoester_all', ['-g', 'all'], alpha_ketoesters)
check('alpha_ketoester_disabled', ['-g', 'all', '-g', '-alpha-ketoester'],
      alpha_ketoesters)
check('alpha_ketoester_enol_preference', ['-g', 'all', '-g', 'ENOL'],
      alpha_ketoesters)

alpha_ketoacids = [
    'O=C(O)C(O)=CC(=O)N', 'O=C(O)C(=O)CC(=O)N',
    'O=C(O)C(O)=CC(=O)Nc1ccccc1', 'O=C(O)C(=O)CC(=O)Nc1ccccc1',
    'O=C(O)C(O)=CC(=O)N(C)C', 'O=C(O)C(=O)CC(=O)N(C)C',
    'O=C(O)C(O)=CC(=O)N1CCCC1', 'O=C(O)C(=O)CC(=O)N1CCCC1',

    'O=C(O)C(O)=Cc1ccccc1', 'O=C(O)C(=O)Cc1ccccc1',
    'O=C(O)C(O)=CC', 'O=C(O)C(O)=C', 'O=C(O)C(O)=Cc1ccncc1',
    '[OH]C(=O)[C]([OH])=[CH]C', '[18OH]C(=O)C(O)=CC',
    'OC(=O)C([18OH])=CC', 'OC(=O)C(O)=[13CH]C',
    'OC(=O)C(O)=C[C@H](F)C', 'O=C([O-])C(O)=CC',
    'O=C(O)C([O-])=CC', 'O=C([O-])C([O-])=CC',
    'O=C(O)/C(O)=C/C', 'O=C(O)C(O[2H])=CC',
    'O=C(O)C(O)=CC(=O)C', 'O=C(O)C(O)=C(C)C',
    'O=C(O)C(O)=CC.O=C(O)C(O)=Cc1ccccc1',
]
check('alpha_ketoacid', ['-g', 'alpha-ketoacid'], alpha_ketoacids)
check('alpha_ketoacid_all', ['-g', 'all'], alpha_ketoacids)
check('alpha_ketoacid_disabled', ['-g', 'all', '-g', '-alpha-ketoacid'],
      alpha_ketoacids)
check('alpha_ketoacid_enol_preference', ['-g', 'all', '-g', 'ENOL'],
      alpha_ketoacids)

biguanides = [
    'CN(C)C(=N)NC(=N)N', 'CN(C)C(N)=NC(=N)N',
    'CN(C)C(=N)N=C(N)N', 'NC(=N)NC(=N)N', 'NC(N)=NC(=N)N',
    'O1CCN(C(=N)NC(=N)N)CC1', 'O1CCN(C(=NC(=N)N)N)CC1',
    'CCNC(=N)NC(=N)N', 'C[C@H](F)NC(=N)NC(=N)N',
    'CN(C)[13C](=N)NC(=N)N', 'CN(C)C(=N)[15NH]C(=N)N',
    'CN(C)C(=N)NC(=[15NH])N', 'CN(C)C(=[NH2+])NC(=N)N',
    'CN(C)C(=N)[NH2+]C(=N)N', 'CN(C)C(N)=[N+]=C(N)N',
    'CN(C)C(=N)N([2H])C(=N)N', 'CCNC(=N)NC(=N)NCC',
    'CN(C)C(=N)NC(=N)N.CN(C)C(=N)NC(=N)N',
]
check('biguanide', ['-g', 'biguanide'], biguanides)
check('biguanide_all', ['-g', 'all'], biguanides)
# Disabled behaviour is covered by the same structural consistency check.
check('biguanide_disabled', ['-g', 'all', '-g', '-biguanide'], biguanides)

fused_lactims = [
    'OC1=NC(O)=CC2=NN(C(N)=C12)c1ccccc1',
    'O=C1NC(O)=CC2=NN(C(N)=C12)c1ccccc1',
    'O=C1NC(=O)C=C2NN(C(N)=C12)c1ccccc1',
    '[18OH]C1=NC(O)=CC2=NN(C(N)=C12)c1ccccc1',
    'OC1=NC([18OH])=CC2=NN(C(N)=C12)c1ccccc1',
    'SC1=NC(S)=CC2=NN(C(N)=C12)c1ccccc1',
    'SC1=NC(O)=CC2=NN(C(N)=C12)c1ccccc1',
    'OC1=NC(S)=CC2=NN(C(N)=C12)c1ccccc1',
    'OC1=NC(O)=CC2=NN(C(N)=C12)[C@H](F)C',
    'OC1=NC(O)=CC2=NN(C(N)=C12)c1ccccc1.CC(O)=NC',
    'OC1=NC(O)=CC2=NN(C(N)=C12)c1ccccc1.'
    'OC1=NC(O)=CC2=NN(C(N)=C12)c1ccccc1',
]
check('fused_lactim', ['-g', 'enol-fused', '-g', 'ltltr'], fused_lactims)
check('fused_lactim_all', ['-g', 'all'], fused_lactims)
check('fused_lactim_disabled', ['-g', 'all', '-g', '-ltltr'], fused_lactims)
check('fused_enol_disabled', ['-g', 'all', '-g', '-enol-fused'], fused_lactims)

guanidines = [
    '[N+]1(=C(N)N)CC2=C(CC1)C=CC=C2',
    '[N+]1(=C(N)N)CCCCC1', 'N=C(N)N1CCCCC1',
    'C[N+](C)=C(N)N', 'C[NH+]=C(N)N',
    '[15N+]1(=C(N)N)CCCCC1', '[N+]1(=[13C](N)N)CCCCC1',
    '[2H][N+](C)=C(N)N', '[N+]1(=C([NH2])[NH2])CCCCC1',
    'C[N+]([C@H](F)C)=C(N)N', 'C[NH+]=C(NCC)NCCC',
    'C[N+](C)=C(NC)NC',
    '[N+]1(=C(N)N)CCCCC1.[N+]1(=C(N)N)CCCCC1',
    '[N+]1(=C(N)N)CCCCC1.C[NH3+]',
]
check('guanidine_charge', ['-g', 'guan'], guanidines)
check('guanidine_charge_all', ['-g', 'all'], guanidines)

amidoximes = [
    'CC(=N)NO', 'CC(N)=NO', 'ONC(=N)C(N)C',
    'ONC(=N)[C@H](N)C', 'ONC(=N)c1ccccc1',
    '[18OH]NC(=N)C', 'O[15NH]C(=N)C', 'ONC(=[15NH])C',
    'ONC(=N)CCC(=N)NO', 'ONC(=N)CCC(=NO)N',
    'ONC(=N)C[N+](C)(C)C', 'CC(=[NH2+])NO',
    'CC(=N)[NH2+]O', 'CC(=N)N[O-]',
    'ON([2H])C(=N)C', 'ONC(=N[2H])C', '[2H]ONC(=N)C',
    'CC(=N)NO.CC(=N)NO',
]
check('amidoxime', ['-g', 'amidoxime'], amidoximes)
check('amidoxime_all', ['-g', 'all'], amidoximes)
check('amidoxime_disabled', ['-g', 'all', '-g', '-amidoxime'], amidoximes)

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
