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
#
#from lillymol import *
#
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

  tmp = Vector{Int32}()
  for i in ring
    push!(tmp, i)
  end
  ring_atoms = OffsetArray(tmp, 0:5)
  # The para exocyclic -O becomes =O
  println("Setting bond btw $(ring[oh_index]) and $(oh)")
  set_bond_type_between_atoms!(m, ring[oh_index], oh, DOUBLE_BOND)

  # Introduce sequence of alternating single and double bonds into `ring`,
  # beginning at the Nitrogen.
  to_place = [1, 2, 1, 1, 2, 1]
  prev = n_index
  println("change: n_index $(n_index) oh_index $(oh_index)")
  # i is the counter and ndx is the index into `ring`.
  for (i,ndx) in enumerate(range(n_index + 1, n_index + 6))
    println("   i $(i) ndx $(ndx)")
    if ndx > 5
      ndx = ndx % 6
    end
    println(" prev $(prev)")
    println(" ndx $(ndx)")
    println(" prev_atom $(ring_atoms[prev])")
    println(" q $(ring_atoms[ndx])")
    println("  prev $(prev) ndx $(ndx) atoms $(ring_atoms[prev]) and $(ring_atoms[ndx])")
    if to_place[i] == 1
      set_bond_type_between_atoms!(m, ring_atoms[prev], ring_atoms[ndx], SINGLE_BOND)
    else
      set_bond_type_between_atoms!(m, ring_atoms[prev], ring_atoms[ndx], DOUBLE_BOND)
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
    if atomic_number(m, a) == 7
      if n_index >= 0
        return
      end
      n_index = ndx
      continue
    end

    println("not N ndx $(ndx) atom $(a) ring $(ring[ndx]) atomic_number $(atomic_number(m, a))")
    # No other heteroatoms
    if atomic_number(m, a) != 6
      continue
    end

    if ncon(m, a) == 2
      continue
    end

    # Look for exocyclic double bonds and the [OD1]
    println("iterate bonds")
    for bond in m[a]
      o = other(bond, a)
      println("  bond to $(o)")
      if o in ring
        continue
      end

      println("SB $(is_single_bond(bond))")
      if ! is_single_bond(bond)
        println("Not single bond, returning")
        return
      end

      println("Exocyclic single bond to type $(atomic_number(m, o))")
      if atomic_number(m, o) != 8
        continue
      end
      if ncon(m, o) != 1
        continue
      end
      if oh_index >= 0
        return
      end
      oh_index = ndx
      oh = o
      println("Just set oh_index $(oh_index) oh $(oh)")
    end
  end

  println("n_index $(n_index) oh_index $(oh_index) carbon $(ring[oh_index]) oh $(oh)")
  if n_index < 0 || oh_index < 0
    return
  end

  # print(f'{ring} N {n_index} c-OH {oh_index} OH {oh}')
  # The two indices must be 3 apart in `ring`. 
  if (n_index + 3) % 6 == oh_index
  elseif (oh_index + 3) == n_index
  else
    return
  end

  println("set bond btw $(ring[oh_index]) and $(oh)")
  println("O bonded $(are_bonded(m, ring[oh_index], oh))")
  change_to_pyridone(m, ring, n_index, oh_index, oh)
end

function four_pyridone(m)
  compute_aromaticity_if_needed(m)
  println("Processing $(name(m))")

  rings = RingAtoms()
  gather_rings(m, rings)
  aromatic = Vector{Bool}()
  for (ndx, r) in enumerate(rings(m))
  aromatic[ndx] = is_aromatic(r)
  end

  for r in rings
    if length(r) != 6
      continue
    end
    println(is_aromatic(r))
    if ! is_aromatic(r)
      continue
    end
    ring_is_four_pyridone(m, r)
  end
  println("Returning four_pyridone")
end


function main(args)
  """Chemical standardisation expt"""
  println(args)
  if length(args) == 0
    println(stderr, "Insufficient arguments")
    return 1
  end
  reader = LillyMol.MoleculeReader(SMI, args[1])
  m = Molecule()
  while next_molecule(reader, m)
    four_pyridone(m)
    molecules_read(reader) > 100000 && break
  end
end

main(ARGS)
