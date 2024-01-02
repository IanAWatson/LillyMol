using LillyMol

function create_fragments(mol)
  natoms(mol) == 1 && return
  compute_aromaticity_if_needed(mol)
  println(smiles(mol))

  # vectors of atom numbers that define the breakable bonds.
  a1s = Vector{Integer}()
  a2s = Vector{Integer}()
  for bond in bonds(mol)
    nrings(bond) > 0 && continue
    is_single_bond(bond) || continue
    push!(a1s, a1(bond))
    push!(a2s, a2(bond))
  end
  length(a1s) == 0 && return

  for (at1, at2) in zip(a1s, a2s)
    remove_bond_between_atoms!(mol, at1, at2)
    for component in create_components(mol)
      create_fragments(component)
    end
    add_bond!(mol, at1, at2, SINGLE_BOND)
  end
end

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

  for component in create_components(mol)
    println(smiles(component))
  end

  mol = LillyMol.MolFromSmiles("CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C([N-]2)C=CC(=C3)OC")
  create_fragments(mol)
  @time create_fragments(mol)
end


main()
