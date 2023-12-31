using LillyMol

function main()
  mol = LillyMol.MolFromSmiles("O-C-C-C-C-N")
  remove_bond_between_atoms!(mol, 0, 1)
  remove_bond_between_atoms!(mol, 2, 3)
  remove_bond_between_atoms!(mol, 4, 5)

  # For each atom that lost a bond, add a dummy atom isotopically
  # labelled with the atom number.
  for atom_number in [0, 1, 2, 3, 4, 5]
    add_atom!(mol, 0)
    add_bond!(mol, atom_number, natoms(mol) - 1, SINGLE_BOND)
    set_isotope!(mol, natoms(mol) - 1, atom_number)
  end

  components = Components()
  create_components(mol, components)
  for component in components
    println(component)
  end
end

main()
