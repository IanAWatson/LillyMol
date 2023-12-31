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
