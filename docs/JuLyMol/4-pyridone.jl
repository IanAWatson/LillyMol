push!(LOAD_PATH, "/home/ian/LillyMol/src/julia")
using ArgMacros
using OffsetArrays
using LillyMol
#
#"""Demo app to show LillyMol python bindings.
#  Convert 4-pyidol to 4-pyidone forms.
#  This was done in order to explore a new chemical standardisation.
#  This implementation was done in python, and when it was found to
#  be doing what we wanted, it was translated to C++.
#
#  Run time for processing the corporate collection with this file
#  was about 6.5 minutes. Doing the same thing with the c++ version
#  was about 1 minute 15 seconds.
#
#  For some tasks this trade-off may be attractive.
#"""

function change_to_pyridone(m, ring, n_index, oh_index, oh)
  """Convert `ring` to 4-pyridone form
    Args:
      m: molecule
      r: 6 membered aromatic ring
      n_index: index in r of [nD2]
      oh_index: index in r of c-[OH]
      oh: atom number of [OH]

    Note that this only works because python has a copy of the rings.
    Once a bond is changed, rings are destroyed. But safe here because
    the ring is not changing.
  """

  # The para exocyclic -O becomes =O
  # println("Setting bond btw $(ring[oh_index]) and $(oh)")
  set_bond_type_between_atoms!(m, ring[oh_index], oh, DOUBLE_BOND)

  # Introduce sequence of alternating single and double bonds into `ring`,
  # beginning at the Nitrogen.
  to_place = [1, 2, 1, 1, 2, 1]
  prev = n_index
  # println("change: n_index $(n_index) oh_index $(oh_index)")
  # i is the counter and ndx is the index into `ring`.
  for (i,ndx) in enumerate(range(n_index + 1, n_index + 6))
    # println("   i $(i) ndx $(ndx)")
    if ndx > 5
      ndx = ndx % 6
    end
    # println(" prev $(prev)")
#     println(" ndx $(ndx)")
 #    println(" prev_atom $(ring[prev])")
  #   println(" q $(ring[ndx])")
   #  println("  prev $(prev) ndx $(ndx) atoms $(ring[prev]) and $(ring[ndx])")
    if to_place[i] == 1
      set_bond_type_between_atoms!(m, ring[prev], ring[ndx], SINGLE_BOND)
    else
      set_bond_type_between_atoms!(m, ring[prev], ring[ndx], DOUBLE_BOND)
    end
    prev = ndx
  end

  println(smiles(m), ' ', name(m))
  true
end

function ring_is_four_pyridone(m, ring)
  """ See if `ring` is a 4-pyridone"""
  # The indices of things in the ring, or the OH
  n_index = -1
  oh_index = -1
  oh = -1;

# Look for [n] and OH in the ring
  for (ndx,a) in enumerate(ring)
    ndx = ndx - 1
    # println("ndx $(ndx) atom $(a)")
    if atomic_number(m, a) == 7
      n_index >= 0 && return

      n_index = ndx
      continue
    end

    # println("not N ndx $(ndx) atom $(a) ring $(ring[ndx]) atomic_number $(atomic_number(m, a))")
    # No other heteroatoms
    atomic_number(m, a) == 6 || continue

    ncon(m, a) == 2 && continue

    # Look for exocyclic double bonds and the [OD1]
    for bond in m[a]
      o = other(bond, a)
      o in ring && continue

      # println("SB $(is_single_bond(bond))")
      is_single_bond(bond) || return

      # println("Exocyclic single bond to type $(atomic_number(m, o))")
      atomic_number(m, o) == 8 || continue
      ncon(m, o) == 1 || continue

      oh_index >= 0 && return

      oh_index = ndx
      oh = o
      # println("Just set oh_index $(oh_index) oh $(oh)")
    end
  end
  if n_index < 0 || oh_index < 0
    return
  end

  # println("n_index $(n_index)")
  # println("oh_index $(oh_index)")
  # println("oh $(oh)")
  # println("N atom $(ring[n_index])")
  # println("carbon $(ring[oh_index])")
  # println("n_index $(n_index) oh_index $(oh_index) carbon $(ring[oh_index]) oh $(oh)")

  # print(f'{ring} N {n_index} c-OH {oh_index} OH {oh}')
  # The two indices must be 3 apart in `ring`. 
  if (n_index + 3) % 6 == oh_index
  elseif (oh_index + 3) == n_index
  else
    return
  end

 #  println("set bond btw $(ring[oh_index]) and $(oh)")
  # println("O bonded $(are_bonded(m, ring[oh_index], oh))")
  change_to_pyridone(m, ring, n_index, oh_index, oh)
end

function four_pyridone(m)
  compute_aromaticity_if_needed(m)
  # println("Processing $(name(m))")

  save_rings = RingInformation()
  gather_rings(m, save_rings, RIP_AROMATIC)
  aromatic = Vector{Bool}(undef, nrings(m))
  for (ndx, r) in enumerate(rings(m))
    aromatic[ndx] = is_aromatic(r)
  end

  for (ndx, r) in enumerate(save_rings)
    # 6 membered aromatics only
    length(r) == 6 || continue
    aromatic[ndx] || continue

    ring_is_four_pyridone(m, r)
  end
#   println("Returning four_pyridone")
end


function main(args)
  """Chemical standardisation expt"""
  if length(args) == 0
    println(stderr, "Insufficient arguments")
    return 1
  end
  reader = LillyMol.MoleculeReader(SMI, args[1])
  m = Molecule()
  while next_molecule(reader, m)
    four_pyridone(m)
    # molecules_read(reader) > 100000 && break
  end
end

main(ARGS)