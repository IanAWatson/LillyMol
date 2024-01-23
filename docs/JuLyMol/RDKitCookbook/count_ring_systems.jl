using LillyMol

function get_ring_systems(mol::Molecule, include_spiro=false)
  # Return a list of the atoms in each ring system.
  if include_spiro
    sysid = label_atoms_by_ring_system_including_spiro_fused(mol)
  else
    sysid = label_atoms_by_ring_system(mol)
  end

  result = []
  nsys = maximum(sysid)

  # non ring atoms get ring system id 0 so start at 1.
  for sys_num in 1:nsys
    push!(result, [i - 1 for (i, s) in enumerate(sysid) if s == sys_num])
  end

  return result
end

function main()
  m = LillyMol.MolFromSmiles("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3")
  rsys = get_ring_systems(m)
  for r in rsys
    println(r)
  end

  m = LillyMol.MolFromSmiles("CC1CC21CC2")
  rsys = get_ring_systems(m, true)
  for r in rsys
    println(r)
  end
end

main()
