using DataStructures

using LillyMol

# Return the Hybridization corresponding to `unsaturation`.
function hybridization(unsaturation)::String
  values = Array{String}([ "SP3", "SP2", "SP"])

  return values[unsaturation + 1]
end

function main()
  m = LillyMol.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)CC#C")

  # Accumulate the number of atoms with each hybridization type.
  total = DefaultDict{String, Integer}(0)
  for i in eachindex(m)
    total[hybridization(unsaturation(m, i))] += 1
  end

  for (h, count) in total
    println("$(h) $(count)")
  end

  # By getting the same information from the atoms.
  total = DefaultDict{String, Integer}(0)
  for atom in m
    total[hybridization(unsaturation(atom))] += 1
  end

  for (h, count) in total
    println("$(h) $(count)")
  end

  # By computing directly.
  total = DefaultDict{String, Int}(0)
  for atom in m
    total[hybridization(nbonds(atom) - ncon(atom))] += 1
  end

  for (h, count) in total
    println("$(h) $(count)")
  end
end

main()
