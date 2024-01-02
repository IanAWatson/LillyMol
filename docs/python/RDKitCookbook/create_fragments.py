from lillymol import *

mol = MolFromSmiles("O-C-C-C-C-N")
mol.remove_bond_between_atoms(0, 1)
mol.remove_bond_between_atoms(2, 3)
mol.remove_bond_between_atoms(4, 5)

# For each atom that lost a bond, add a dummy atom isotopically
# labelled with the atom number.
for atom_number in [0, 1, 2, 3, 4, 5]:
  mol.add_atom(0)
  mol.add_bond(atom_number, mol.natoms() - 1, SINGLE_BOND)
  mol.set_isotope(mol.natoms() - 1, atom_number)

for component in mol.create_components():
  print(component)

def create_fragments(mol:Molecule):
  mol.compute_aromaticity_if_needed()
  print(mol.smiles())
  breakable_bonds = []

  for bond in mol.bonds():
    if bond.nrings():
      continue
    if not bond.is_single_bond():
      continue
    breakable_bonds.append(bond)

  for bond in breakable_bonds:
    mol.remove_bond_between_atoms(bond.a1(), bond.a2())
    for component in mol.create_components():
      create_fragments(component);
    mol.add_bond(bond.a1(), bond.a2(), SINGLE_BOND)

mol = MolFromSmiles("CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C([N-]2)C=CC(=C3)OC")
create_fragments(mol)


