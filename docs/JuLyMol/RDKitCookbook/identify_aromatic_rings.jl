using LillyMol

function main()
  mol = LillyMol.MolFromSmiles("c1cccc2c1CCCC2")
  compute_aromaticity_if_needed(mol)
  for ring in rings(mol)
    println(ring)
    println(is_aromatic(ring))
  end
end

main()
