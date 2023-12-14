module LillyMol
  using CxxWrap

  # abstract type AbstractSetOfAtoms <: AbstractVector{Int32} end

  import Base: getindex, iterate, in, length, size, push!, show

  # @wrapmodule(joinpath("bazel-bin/julia/","lillymol_julia.so"))
  @wrapmodule(() -> joinpath("bazel-bin/julia/","lillymol_julia.so"))
#  @wrapmodule(libfoo_jll.get_libfoo_path

  function __init__()
    @initcxx
  end
  export World
  export greet

  export FileType
  export BondType, SINGLE, DOUBLE, TRIPLE, AROMATIC

  export Molecule, SetOfAtoms, Atom, Bond, ChemicalStandardisation, BondList, Mol2Graph, ChiralCentre
  export SetOfRings

  # Now done in C++
  # getindex(m::Molecule, a::Int)=atom(m, a)
  # getindex(a::Atom, b::Int)=atom.item(b)
  iterate(m::Molecule, state=0) = (state >= natoms(m) ? nothing : (m[state], state + 1))
  iterate(a::Atom, state=0) = (state >= ncon(a) ? nothing : (a[state], state + 1))
  iterate(b::Bond, state=1) = (state == 1 ? (b.a1(), 2) : state == 2 ? (b.a2(), 2) : nothing)
  iterate(b::BondList, state=0) = (state >= bonds_in_set(b) ? nothing : (b[state], state + 1))
  iterate(r::Ring, state=0) = (state >= length(r) ? nothing : (r[state], state + 1))
  iterate(s::SetOfAtoms, state=0) = (state >= length(s) ? nothing : (s[state], state + 1))
  iterate(r::SetOfRings, state=0) = (state >= length(r) ? nothing : (r[state], state + 1))
  in(z::Int, m::Molecule) = (natoms(m, z) > 0)
  in(atom::Int, a::Atom) = involves(a, atom)
  length(m::Molecule) = natoms(m)
  length(r::Ring) = atoms_in_ring(r)
  length(s::SetOfRings) = rings_in_set(s)
  length(b::BondList) = bonds_in_set(b)
  # length(r::Ring) = size(r)
  size(m::Molecule) = natoms(m)
  size(r::Ring) = (atoms_in_ring(r),)
  size(s::SetOfAtoms) = (length(s),)
  size(s::SetOfAtomsAllocated) = (length(s),)
  size(s::SetOfRings) = (rings_in_set(s),)
  size(s::SetOfAtomsAllocated) = (length(s),)
  size(b::BondList) = (length(b),)
  export getindex
  export iterate
  export length
  export in
  export atoms_in_ring, contains
  export natoms, smiles, unique_smiles, nrings, build_from_smiles, is_aromatic, set_name!, name
  export atomic_number, molecular_formula, nedges, is_ring_atom, fused_system_size, fused_system_identifier
  export rings_with_fused_system_identifier, in_same_ring, in_same_aromatic_ring, in_same_ring_system
  export ring_membership, rings_containing_both, is_part_of_fused_ring_system, ring, ring_containing_atom
  export label_atoms_by_ring_system, label_atoms_by_ring_system_including_spiro_fused
  export nrings_including_non_sssr_rings, non_sssr_rings, non_sssr_ring, is_spiro_fused, is_halogen
  export maximum_connectivity, connections, isotopically_labelled_smiles, is_aromatic, atom, formal_charge
  export isotope, number_formally_charged_atoms, net_formal_charge, bond, bond_between_atoms
  export compute_aromaticity_if_needed, bond_between_atoms, number_symmetry_classes, symmetry_class, symmetry_equivalents
  export symmetry_classes, attached_heteroatom_count, add_bond, are_bonded, bond_between_atoms
  export sssr_rings, rings_in_set, bond_list, bonds_in_set, add!, remove_atom!, remove_atoms!, delete_fragment!
  export remove_fragment_containing_atom, remove_all, atomic_symbol, remove_all_non_natural_elements
  export remove_explicit_hydrogens, valence_ok, remove_bonds_to_atom, remove_bond, remove_bond_between_atoms
  export remove_all_bonds, molecular_weight, molecular_weight_count_isotopes, molecular_weight_ignore_isotopes
  export bond_length, bond_angle, dihedral_angle, highest_coordinate_dimensionality, exact_mass
  export translate, discern_chirality_from_3d_structure
  export centre, top_front, top_back, left_down, right_down
  export chiral_centres, chiral_centre_at_atom, chiral_centre_in_molecule_not_indexed_by_atom_number
  export remove_chiral_centre_at_atom, remove_all_chiral_centres, invert_chirality_on_atom
  export number_fragments, fragment_membership, atoms_in_fragment, get_atoms_in_fragment, largest_fragment
  export identify_spinach, rings_in_fragment, create_components, create_subset
  export reduce_to_largest_fragment, reduce_to_largest_organic_fragment, reduce_to_largest_fragment_carefully
  export organic_only, contains_non_periodic_table_elements
  export longest_path, atoms_between, bonds_between
  export implicit_hydrogens, explicit_hydrogens, hcount, move_hydrogens_to_end_of_connection_table
  export make_implicit_hydrogens_explicit, pi_electrons, lone_pair_count, saturated
  export aromatic_atom_count, aromatic_ring_count, unset_unnecessary_implicit_hydrogens_known_values
  export smarts_equivalent_for_atom, smarts
  export atom_map_number, set_atom_map_number, atom_with_atom_map_number, reset_all_atom_map_numbers
  export set_include_atom_map_with_smiles
  export x, y, z, set_x, setx, sety, setz, setxyz
  export ncon, nbonds, involves, other
  export is_single_bond, is_double_bond, is_triple_bond, is_aromatic
  export atomic_number
  export a1, a2
  show(io::IO, s::SetOfAtoms) = print(io, set_of_atoms_show_text(s))
  # For some reason this does not work, but works for set_of_atoms.
  show(io::IO, r::Ring) = print(io, ring_show_text(r))
  show(io::IO, m::Molecule) = print(io, molecule_show_text(m))

  export activate_all, process
  # Not sure why this was necessary. The returned BondList from bond_list(m) did not
  # generate a proper getindex method, and defining this method was suggested
  # as the solution.
  getindex(blist::Union{CxxWrap.CxxWrapCore.ConstCxxRef{<:BondList}, CxxWrap.CxxWrapCore.CxxRef{<:BondList}}, ndx::Int64) =
    internal_get_item(blist, ndx)
end

function boobar()
  return true
  println("boobar begins")
  m = LillyMol.MolFromSmiles("CCC.C1CC1.C1CC1C2CC2")
  nf = number_fragments(m)
  nf == 3 || return false
  frags = [Molecule() for i in 1:nf]
  create_components(m, frags)
  length(frags) == 3 || return false
  println("boobar done")
  true
end
#getindex(m::LillyMol.Molecule, a::Int64)=LillyMol.atom(m, a)

using .LillyMol
# Never a great idea to use random with tests...
using Random
using Test

# Call greet and show the result
@show LillyMol.greet()

w = LillyMol.World()
println(LillyMol.greet(w));
LillyMol.set(w, "hello")
println(LillyMol.greet(w))

# Mimic the same functionality from GoogleTest
# Implemented after I had finished this, use for new cases.
function unordered_elements_are_array(v1::Vector, v2::Vector)::Bool
  size(v1) != size(v2) && return false
  return sort(v1) == sort(v2)
end

function test_set_of_atoms()::Bool
  s = SetOfAtoms()
  for i in 1:6
    push!(s, i)
  end
  for i in 1:6
    i in s || return false
  end
  for i in 1:6
    s[i - 1] = 6 - i
  end
  collect(s) == [5, 4, 3, 2, 1, 0] || return false

  true
end

function test_empty_molecule()::Bool
  m = Molecule()
  natoms(m) == 0 || return false
  length(name(m)) == 0 || return false
  return true
end

function test_build_from_smiles()::Bool
  m = Molecule()
  build_from_smiles(m, "c1ccccc1OC foo")
end

function test_set_name()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  set_name!(m, "foo")
  name(m) == "foo"
end

function test_natoms()::Bool
  m = Molecule()
  for i in 1:20
    smiles = "C" ^ i
    build_from_smiles(m, smiles) || return false
    natoms(m) == i || return false
  end
  true
end

function test_natoms_atomic_number()::Bool
  m = Molecule()
  for i in 1:10
    smiles = "CN" ^ i
    build_from_smiles(m, smiles) || return false
    natoms(m, 6) == i || return false
    natoms(m, 7) == i || return false
  end
  true
end

function test_natoms_element()::Bool
  m = Molecule()
  for i in 1:10
    smi = "CN" ^ i
    build_from_smiles(m, smi) || return false
    natoms(m, "C") == i || return false
    natoms(m, "N") == i || return false
  end
  true
end

function test_length_molecule()::Bool
  m = Molecule()
  for i = 1:10
    smi = "CN"^i
    build_from_smiles(m, smi) || return false
    length(m) == (2 * i) || return false
  end
  true
end

function test_atomic_number()::Bool
  atomic_numbers = [9, 6, 7, 8, 92]
  m = Molecule()
  build_from_smiles(m, "FCNO.[U]")
  for (ndx, atom) in enumerate(m)
    atomic_number(m, ndx - 1) == atomic_numbers[ndx] || return false
    atomic_number(atom) == atomic_numbers[ndx] || return false
  end
  true
end

function test_nedges()::Bool
  mol_and_edges = Dict{String, Int}(""=>0, "[Fe]"=>0, "II"=>1, "III"=>2, "I1II1"=>3)
  for (smi,edges) in mol_and_edges
    m = Molecule()
    build_from_smiles(m, smi) || return false
    nedges(m) == edges || "$(edges) failed" == ""
  end
  true
end

function test_atomic_number_in():Bool
  m = Molecule()
  for i in 1:10
    smi = "CNOC(F)(F)" ^ i
    build_from_smiles(m, smi) || return false
    for z in [6, 7, 8, 9]
      z in m || return false
    end
  end
  true
end

function test_molecular_formula()::Bool
  smi_to_mf = Dict{String, String}("C"=>"CH4", "CO"=>"CH4O",
        "C.C"=>"CH4.CH4", "c1ccccc1"=>"C6H6")
  for (smi, mf) in smi_to_mf
    m = Molecule()
    build_from_smiles(m, smi) || return false
    # println("Expct $(mf) get $(molecular_formula(m))")
    mf == molecular_formula(m) || return false
  end
  true
end

function test_valence_ok()::Bool
  m = Molecule()
  build_from_smiles(m, "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F") || return false
  valence_ok(m) || return false
  build_from_smiles(m, "[CH2]-C") || return false
  valence_ok(m) && return false
  build_from_smiles(m, "CC(C)(C)(C)(C)(C)(C)(C)(C)(C)C") || return false
  valence_ok(m) && return false
  true
end

function test_nrings()::Bool
  smiles_and_nrings = Dict{String, Int}("C"=>0, "CC"=>0, "CCC"=>0,
        "C1CC1"=>1, "C12CC1C2"=>2, "C1CC1C1CC1"=>2)
  m = Molecule()
  for (smi, nr) in smiles_and_nrings
    build_from_smiles(m, smi) || return false
    nrings(m) == nr || return false
  end
  true
end

function test_nrings_atom()::Bool
  smiles_and_nrings = Dict{String, Array{Int}}("C"=>[0], "CC"=>[0,0], "CCC"=>[0,0,0],
        "C1CC1"=>[1,1,1], "C12CC1C2"=>[2,1,2,1], "C1CC1C1CC1"=>[1,1,1,1,1,1])
  m = Molecule()
  for (smi, rings) in smiles_and_nrings
    build_from_smiles(m, smi) || return false
    for (atom, nr) in enumerate(rings)
      atom = atom - 1
      nrings(m, atom) == nr || return false
    end
  end
  true
end

function test_nrings_size()::Bool
  m = Molecule()
  for rsize in 3:10
    smi = "C1" * "C"^(rsize-1) * "1"
    build_from_smiles(m, smi) || return false
    for i in 1:natoms(m)
      nrings(m, i - 1, rsize) == 1 || return false
    end
  end
  build_from_smiles(m, "C12CC1C2") || return false
  nrings(m, 0, 3) == 2 || return false
  nrings(m, 1, 3) == 1 || return false
  nrings(m, 2, 3) == 2 || return false
  nrings(m, 3, 3) == 1 || return false
  true
end

function test_is_ring_atom()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  is_ring_atom(m, 0) == false || return false
  build_from_smiles(m, "C1CC1C")
  is_ring_atom(m, 0) == true || return false
  is_ring_atom(m, 1) == true || return false
  is_ring_atom(m, 2) == true || return false
  is_ring_atom(m, 3) == false || return false

  true
end

function test_ring_bond_count()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  ring_bond_count(m, 1) == 0 || return false
  build_from_smiles(m, "C1CC1C")
  ring_bond_count(m, 0) == 1 || return false
  ring_bond_count(m, 1) == 1 || return false
  ring_bond_count(m, 2) == 1 || return false
  ring_bond_count(m, 3) == 0 || return false

  build_from_smiles(m, "C12CC2C1") || return false
  ring_bond_count(m, 0) == 2 || return false
  ring_bond_count(m, 1) == 1 || return false
  ring_bond_count(m, 2) == 2 || return false
  ring_bond_count(m, 3) == 1 || return false

  true
end

function test_fused_system_size()::Bool
  m = Molecule();
  build_from_smiles(m, "C") || return false
  fused_system_size(m, 0) == 0 || return false
  build_from_smiles(m, "C1CC1C") || return false
  fused_system_size(m, 0) == 1 || return false
  fused_system_size(m, 3) == 0 || return false

  build_from_smiles(m, "C12CC2C1") || return false
  for i in 1:natoms(m)
    fused_system_size(m, i - 1) == 2 || return false
  end
  true
end

function test_fused_system_identifier()::Bool
  m = Molecule()
  build_from_smiles(m, "CC") || return false
  fused_system_identifier(m, 0) == fused_system_identifier(m, 1) || return false
  build_from_smiles(m, "C1CC1") || return false
  fused_system_identifier(m, 0) == fused_system_identifier(m, 1) || return false
  fused_system_identifier(m, 1) == fused_system_identifier(m, 2) || return false
  build_from_smiles(m, "C12CC2C1") || return false
  for i in 2:natoms(m)
    fused_system_identifier(m, 0) == fused_system_identifier(m, i - 1) || return false
  end
  true
end

function test_rings_with_fused_system_identifier()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1")
  fsid = fused_system_identifier(m, 0);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  build_from_smiles(m, "C1CC1C1CC1")
  fsid = fused_system_identifier(m, 0);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  fsid = fused_system_identifier(m, 3);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  build_from_smiles(m, "C12CC2C1") || return false
  fsid = fused_system_identifier(m, 3);
  rings_with_fused_system_identifier(m, fsid) == 2 || return false
  true
end

function test_in_same_ring()::Bool
  m = Molecule();
  build_from_smiles(m, "CC") || return false
  in_same_ring(m, 0, 1) && return false
  build_from_smiles(m, "C1CC1C") || return false
  in_same_ring(m, 0, 1) || return false
  in_same_ring(m, 1, 2) || return false
  in_same_ring(m, 0, 2) || return false
  ! in_same_ring(m, 0, 3) || return false
  ! in_same_ring(m, 1, 3) || return false
  ! in_same_ring(m, 2, 3) || return false

  build_from_smiles(m, "C12CC2C1") || return false
  in_same_ring(m, 0, 1) || return false
  in_same_ring(m, 0, 2) || return false
  in_same_ring(m, 0, 3) || return false
  in_same_ring(m, 1, 2) || return false
  ! in_same_ring(m, 1, 3) || return false
  in_same_ring(m, 2, 3) || return false
end

function test_in_same_aromatic_ring()::Bool
  m = Molecule()
  build_from_smiles(m, "Oc1cccnc1") || return false
  for i in 1:6
    in_same_aromatic_ring(m, 0, i - 1) && return false
  end
  for i in 1:6
    for j in i:6
      in_same_aromatic_ring(m, i, j) || return false
    end
  end
  build_from_smiles(m, "c12cncc2cccc1") || return false
  in_same_aromatic_ring(m, 1, 2) || return false
  in_same_aromatic_ring(m, 1, 3) || return false
  in_same_aromatic_ring(m, 1, 4) || return false
  ! in_same_aromatic_ring(m, 1, 5) || return false
  ! in_same_aromatic_ring(m, 5, 1) || return false
  ! in_same_aromatic_ring(m, 1, 6) || return false
  ! in_same_aromatic_ring(m, 6, 1) || return false
  ! in_same_aromatic_ring(m, 1, 7) || return false
  ! in_same_aromatic_ring(m, 1, 8) || return false
end

function test_in_same_ring_system()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1CC1CC1") || return false
  in_same_ring_system(m, 0, 1) || return false
  in_same_ring_system(m, 0, 2) || return false
  in_same_ring_system(m, 1, 2) || return false
  in_same_ring_system(m, 4, 5) || return false
  in_same_ring_system(m, 4, 6) || return false
  in_same_ring_system(m, 5, 6) || return false
  ! in_same_ring_system(m, 2, 3) || return false
  ! in_same_ring_system(m, 3, 4) || return false
  ! in_same_ring_system(m, 2, 4) || return false
  ! in_same_ring_system(m, 1, 5) || return false
  ! in_same_ring_system(m, 0, 6) || return false
end

function test_ring_membership()::Bool
  m = Molecule();
  build_from_smiles(m, "CCC") || return false
  r = ring_membership(m)
  sum(r) == 0 || return false
  build_from_smiles(m, "C1CC1") || return false
  r = ring_membership(m)
  # https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal#:~:text=The%20shortest%20way%20I%20can,%3D%3D%20arr)%20.
  # all(y->y==r[1], r) || return false
  all(y->y==1, r) || return false
  build_from_smiles(m, "C12CC2C1")
  r = ring_membership(m)
  [2, 1, 2, 1] == r || return false
  true
end

function test_rings_containing_both()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C") || return false
  rings_containing_both(m, 0, 1) == 1 || return false
  rings_containing_both(m, 0, 2) == 1 || return false
  rings_containing_both(m, 1, 2) == 1 || return false
  rings_containing_both(m, 0, 3) == 0 || return false
  rings_containing_both(m, 1, 3) == 0 || return false
  rings_containing_both(m, 2, 3) == 0 || return false
  build_from_smiles(m, "C12CC2C1C") || return false
  rings_containing_both(m, 0, 1) == 1 || return false
  rings_containing_both(m, 0, 2) == 2 || return false
  rings_containing_both(m, 0, 3) == 1 || return false
  rings_containing_both(m, 1, 2) == 1 || return false
  rings_containing_both(m, 2, 3) == 1 || return false
  rings_containing_both(m, 3, 4) == 0 || return false
  true
end

function test_is_part_of_fused_ring_system()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return false
  ! is_part_of_fused_ring_system(m, 0) || return false
  build_from_smiles(m, "C1CC1") || return false
  ! is_part_of_fused_ring_system(m, 0) || return false
  ! is_part_of_fused_ring_system(m, 1) || return false
  build_from_smiles(m, "C12CC2C1CC") || return false
  is_part_of_fused_ring_system(m, 0) || return false
  is_part_of_fused_ring_system(m, 1) || return false
  is_part_of_fused_ring_system(m, 2) || return false
  is_part_of_fused_ring_system(m, 3) || return false
  ! is_part_of_fused_ring_system(m, 4) || return false
end

function test_ring()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C") || return false
  r = ring(m, 0)
  l = length(r)
  length(r) == 3 || return false
  for i in 0:2
    i in r || return false
  end
  atoms = collect(r)
  for i in 0:2
    i in atoms || return false
  end
  3 in atoms && return false

  build_from_smiles(m, "C1CC1CC1CC1") || return false
  for i in 1:nrings(m)
    r = ring(m, i - 1)
    length(r) == 3 || return false
  end
  all([x in ring(m, 0) for x in 0:2]) || return false
  all([x in ring(m, 1) for x in 4:6]) || return false
  true
end

function same_atoms(s1, s2)::Bool
  length(s1) == length(s2) || return false
  set1 = Set{Int}()
  for i in s1
    push!(set1, i)
  end
  for i in s2
    i in set1 || return false
  end
  true
end

function test_rings()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1C2CCC2C3CCCC3") || return false
  expected = [[0, 1, 2], [3, 4, 5, 6], [7, 8, 9, 10, 11]]
  for (i, r) in enumerate(sssr_rings(m))
    same_atoms(expected[i], r) || return false
  end
  true
end

function test_ring_containing_atom()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C2CCC2C3CCCC3") || return false
  for i in 0:2
    r = ring_containing_atom(m, i)
    length(r) == 3 || return false
  end
  for i in 3:6
    r = ring_containing_atom(m, i)
    length(r) == 4 || return false
  end
  for i in 7:(natoms(m) - 1)
    r = ring_containing_atom(m, i)
    length(r) == 5 || return false
  end
  true
end

# Fix for zero arrays
function test_label_atoms_by_ring_system()::Bool
  m = LillyMol.MolFromSmiles("C1CC1C2CCC2C3CCCCC3C")
  rsys = label_atoms_by_ring_system(m)
  for i in 1:3
    rsys[i] == 1 || return false
  end
  for i in 4:7
    rsys[i] == 2 || return false
  end
  for i in 8:13
    rsys[i] == 3 || return false
  end
  rsys[14] == 0 || return false

  true
end

# Fix for zero arrays
function test_label_atoms_by_ring_system_including_spiro_fused()::Bool
  m = LillyMol.MolFromSmiles("C1CC12CC2")
  rsys = label_atoms_by_ring_system_including_spiro_fused(m)
  for i in 1:natoms(m)
    rsys[i] == 1 || return false
  end
  true
end

function test_nrings_including_non_sssr_rings()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  for i in 0:(natoms(m) - 1)
    nrings_including_non_sssr_rings(m, i) == 3 || return false
  end
  true
end

function test_non_sssr_rings()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  non_sssr_rings(m) == 1 || return false
  true
end

function test_non_sssr_ring()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  # 1 == non_sssr_rings(m) || return false
  r = non_sssr_ring(m, 0)
  length(r) == 4 || return false
  true
end

function test_is_spiro_fused()::Bool
  m = LillyMol.MolFromSmiles("C1CC12CC2")
  is_spiro_fused(m, 0) && return false
  is_spiro_fused(m, 1) && return false
  is_spiro_fused(m, 2) || return false
  is_spiro_fused(m, 3) && return false
  is_spiro_fused(m, 4) && return false
  true
end

function test_is_halogen()::Bool
  m = LillyMol.MolFromSmiles("FC(Cl)(Br)CI")
  is_halogen(m, 0) || return false
  is_halogen(m, 1) && return false
  is_halogen(m, 2) || return false
  is_halogen(m, 3) || return false
  is_halogen(m, 4) && return false
  is_halogen(m, 5) || return false
  true
end

function test_ncon_molecule()::Bool
  m = LillyMol.MolFromSmiles("CCC(C)C(C)(C)C")
  ncon(m, 0) == 1 || return false
  ncon(m, 1) == 2 || return false
  ncon(m, 2) == 3 || return false
  ncon(m, 4) == 4 || return false
  true
end

function test_nbonds_molecule()::Bool
  m = LillyMol.MolFromSmiles("CCC=CC#C")
  nbonds(m, 0) == 1 || return false
  nbonds(m, 1) == 2 || return false
  nbonds(m, 2) == 3 || return false
  nbonds(m, 4) == 4 || return false
  true
end

function test_maximum_connectivity()::Bool
  m = Molecule()
  maximum_connectivity(m) == 0 || return false
  build_from_smiles(m, "C") || return false
  maximum_connectivity(m) == 0 || return false
  build_from_smiles(m, "CC") || return false
  maximum_connectivity(m) == 1 || return false
  build_from_smiles(m, "CCC") || return false
  maximum_connectivity(m) == 2 || return false
  build_from_smiles(m, "CC(C)C") || return false
  maximum_connectivity(m) == 3 || return false
  build_from_smiles(m, "CC(C)(C)C") || return false
  maximum_connectivity(m) == 4 || return false
  true
end

function test_connections_molecule()::Bool
  m = Molecule();
  build_from_smiles(m, "C") || return false
  c = connections(m, 0)
  length(c) == 0 || return false
  build_from_smiles(m, "CC") || return false
  c = connections(m, 0)
  length(c) == 1 || return false
  c[1] == 1 || return false
  c = connections(m, 1)
  length(c) == 1 || return false
  c[1] == 0 || return false
  build_from_smiles(m, "CCC") || return false
  c = connections(m, 1)
  (0 in c && 2 in c) || return false
  build_from_smiles(m, "CC(C)(C)C") || return false
  c = connections(m, 1)
  length(c) == 4 || return false
  (0 in c && 2 in c && 3 in c && 4 in c) || return false
  unordered_elements_are_array(collect(c), [0, 2, 3, 4]) || return false
  true
end

function test_isotopically_labelled_smiles()::Bool
  m = Molecule()
  build_from_smiles(m, "CNO")
  s = isotopically_labelled_smiles(m)
  s == "C[1NH][2OH]" || return false
  true
end

function test_is_aromatic()::Bool
  m = LillyMol.MolFromSmiles("C1CCC1")
  for i in natoms(m)
    is_aromatic(m, i - 1) && return false
  end

  build_from_smiles(m, "c1ccccc1") || return false
  for i in natoms(m)
    is_aromatic(m, i - 1) || return false
  end
  true
end

function test_getindex_molecule()::Bool
  m = LillyMol.MolFromSmiles("[1CH3]NOF")
  atomic_number(m[0]) == 6 || return false
  atomic_number(m, 0) == 6 || return false
  atomic_number(m[1]) == 7 || return false
  atomic_number(m, 1) == 7 || return false
  atomic_number(m[2]) == 8 || return false
  atomic_number(m, 2) == 8 || return false
  atomic_number(m[3]) == 9 || return false
  atomic_number(m, 3) == 9 || return false
  atomic_symbol(m[0]) == "C" || return false
  atomic_symbol(m, 0) == "C" || return false
  ncon(m[0]) == 1 || return false
  ncon(m, 0) == 1 || return false
  ncon(m[1]) == 2 || return false
  ncon(m, 1) == 2 || return false
  ncon(m[2]) == 2 || return false
  ncon(m, 2) == 2 || return false
  ncon(m[3]) == 1 || return false
  ncon(m, 3) == 1 || return false
  formal_charge(m[0]) == 0 || return false
  formal_charge(m, 0) == 0 || return false
  isotope(m[0]) == 1 || return false
  isotope(m, 0) == 1 || return false
  isotope(m[1]) == 0 || return false
  isotope(m, 1) == 0 || return false
  for (ndx, atom) in enumerate(m)
    ndx = ndx - 1
    atomic_number(m[ndx]) == atomic_number(atom) || return false
    atomic_symbol(m[ndx]) == atomic_symbol(atom) || return false
    ncon(m[ndx]) == ncon(atom) || return false
    formal_charge(m[ndx]) == formal_charge(atom) || return false
    isotope(m[ndx]) == isotope(atom) || return false
  end
  true
end

function test_number_formally_charged_atoms()::Bool
  m = LillyMol.MolFromSmiles("NC")
  number_formally_charged_atoms(m) == 0 || return false
  build_from_smiles(m, "[NH3+]C")
  number_formally_charged_atoms(m) == 1 || return false
  build_from_smiles(m, "C(=O)[O-]")
  number_formally_charged_atoms(m) == 1 || return false
  true
end

function test_net_formal_charge()::Bool
  m = LillyMol.MolFromSmiles("NC")
  net_formal_charge(m) == 0 || return false
  build_from_smiles(m, "[NH3+]C")
  net_formal_charge(m) == 1 || return false
  build_from_smiles(m, "C(=O)[O-]")
  net_formal_charge(m) == -1 || return false
  true
end

function test_bond_molecule()::Bool
  m = LillyMol.MolFromSmiles("CC=CC#C")
  b = bond(m, 0)
  is_single_bond(b) || return false
  is_double_bond(bond(m, 1)) || return false
  is_triple_bond(bond(m, 3)) || return false
  build_from_smiles(m, "c1ccccc1")
  compute_aromaticity_if_needed(m)
  for i in 1:nedges(m)
    is_aromatic(bond(m, i - 1)) || return false
  end
  true
end

function test_bond_between_atoms()::Bool
  m = LillyMol.MolFromSmiles("CC=C")
  b = bond_between_atoms(m, 0, 1)
  is_single_bond(b) || return false
  b = bond_between_atoms(m, 1, 2)
  is_double_bond(b) || return false
  true
end

function test_number_symmetry_classes()::Bool
  m = LillyMol.MolFromSmiles("C")
  number_symmetry_classes(m) == 1 || return false
  build_from_smiles(m, "CC") || return false
  number_symmetry_classes(m) == 1 || return false
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  number_symmetry_classes(m) == 4 || return false
  build_from_smiles(m, "C1CC1") || return false
  number_symmetry_classes(m) == 1 || return false
  true
end

function test_symmetry_class()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  symmetry_class(m, 0) == symmetry_class(m, 2) || return false
  symmetry_class(m, 0) == symmetry_class(m, 3) || return false
  symmetry_class(m, 5) == symmetry_class(m, 6) || return false
  build_from_smiles(m, "CN")
  symmetry_class(m, 0) == symmetry_class(m, 1) && return false
  true
end

function test_symmetry_equivalents()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  s = symmetry_equivalents(m, 0)
  length(s) == 2 || return false
  (2 in s && 3 in s) || return false
  s = symmetry_equivalents(m, 6)
  length(s) == 1 || return false
  (6 in s) || return false
  true
end

# Fix for zero arrays
function test_symmetry_classes()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  c = symmetry_classes(m)
  c[1] == c[3] || return false
  c[1] == c[4] || return false
  c[6] == c[7] || return false
  true
end

function test_attached_heteroatom_count()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)N") || return false
  attached_heteroatom_count(m, 0) == 0 || return false
  attached_heteroatom_count(m, 1) == 3 || return false
  attached_heteroatom_count(m, 2) == 0 || return false
  attached_heteroatom_count(m, 4) == 1 || return false
  attached_heteroatom_count(m, 5) == 0 || return false
  attached_heteroatom_count(m, 6) == 0 || return false
  true
end

function test_bond_length()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,1,1}}") || return false
  @test bond_length(m, 0, 1) ≈ sqrt(3.0) atol=0.001
  true
end

function test_bond_angle()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,0,0}}C{{1,1,1}}") || return false
  @test bond_angle(m, 0, 1, 2) ≈ (π / 2.0) atol=0.001
  true
end

function test_dihedral_angle()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,0,0}}C{{1,1,1}}") || return false
  @test bond_angle(m, 0, 1, 2) ≈ (π / 2.0) atol=0.001
  true
end

function test_add_bond()::Bool
  m = Molecule()
  build_from_smiles(m, "C.C") || return false
  add_bond(m, 0, 1, SINGLE)
  smiles(m) == "CC" || return false
  build_from_smiles(m, "C.C") || return false
  add_bond(m, 0, 1, DOUBLE)
  smiles(m) == "C=C" || return false
  build_from_smiles(m, "C.C") || return false
  add_bond(m, 0, 1, TRIPLE)
  smiles(m) == "C#C" || return false
  true
end

function test_are_bonded()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC")
  are_bonded(m, 0, 1) && return false
  are_bonded(m, 1, 2) || return false
  true
end

function test_bond_between_atoms()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return false
  b = bond_between_atoms(m, 0, 1)
  ((a1(b) == 0 && a2(b) == 1) || (a1(b) == 1 && a2(b) == 0)) || return false
  true
end

function test_bond_list()::Bool
  m = Molecule()
  build_from_smiles(m, "CC=CC#C") || return false
  for (bnum, bond) in enumerate(bond_list(m))
    if bnum == 1
      (a1(bond) == 0 && a2(bond) == 1) || return false
    elseif bnum == 2
      (a1(bond) == 1 && a2(bond) == 2) || return false
    elseif bnum == 3
      (a1(bond) == 2 && a2(bond) == 3) || return false
    elseif bnum == 4
      (a1(bond) == 3 && a2(bond) == 4) || return false
    end
  end
  for (bnum, bond) in enumerate(bond_list(m))
    if bnum == 1
      is_single_bond(bond) || return false
    elseif bnum == 2
      is_double_bond(bond) || return false
    elseif bnum == 3
      is_single_bond(bond) || return false
    elseif bnum == 4
      is_triple_bond(bond) || return false
    end
  end
  true
end

function test_add_molecule()::Bool
  m1 = LillyMol.MolFromSmiles("C")
  m2 = LillyMol.MolFromSmiles("N")
  add!(m1, m2)
  smiles(m1) == "C.N" || return false
  add!(m1, m2)
  smiles(m1) == "C.N.N" || return false
  true
end

function test_remove_atom()::Bool
  m = LillyMol.MolFromSmiles("CN")
  remove_atom!(m, 0)
  smiles(m) == "N" || return false
  true
end

function test_remove_atoms()::Bool
  m = LillyMol.MolFromSmiles("CNCOF")
  s = SetOfAtoms()
  add!(s, 1)
  add!(s, 3)
  remove_atoms!(m, s)
  smiles(m) == "C.C.F" || return false
  true
end

function test_remove_atoms_vector()::Bool
  m = LillyMol.MolFromSmiles("CNCOF")
  to_remove = [0, 1, 0, 1, 0]
  remove_atoms!(m, to_remove) == 2 || return false
  smiles(m) == "C.C.F" || return false
  true
end

function test_delete_fragment()::Bool
  m = LillyMol.MolFromSmiles("CCC.NNN.O")
  delete_fragment!(m, 1)
  smiles(m) == "CCC.O" || return false
  true
end

function test_remove_fragment_containing_atom()::Bool
  m = LillyMol.MolFromSmiles("C.N.O.S")
  remove_fragment_containing_atom(m, 1)
  smiles(m) == "C.O.S" || return false
  true
end

function test_remove_all()::Bool
  m = LillyMol.MolFromSmiles("OC(=O)C")
  remove_all(m, 8) == 2 || return false
  smiles(m) == "CC" || return false
  true
end

function test_set_auto_create_new_elements()::Bool
  LillyMol.set_auto_create_new_elements(1)
  m = LillyMol.MolFromSmiles("[Th][Eq][U]IC[K][Br]O[W]NFO[Xj][Um]PSO[Ve][Rt][La][Zy][D]O[G]")
  length(m) == 24 || return false
  true
end

function test_set_atomic_symbols_can_have_arbitrary_length()::Bool
  LillyMol.set_atomic_symbols_can_have_arbitrary_length(1)
  m = LillyMol.MolFromSmiles("[Ala][Gly][Ser][Hello]")
  natoms(m) == 4 || return false
  atomic_symbol(m, 0) == "Ala" || return false
  atomic_number(m, 1) < 0 || return false
  LillyMol.set_atomic_symbols_can_have_arbitrary_length(0)
  true
end

function test_remove_all_non_natural_elements()::Bool
  LillyMol.set_auto_create_new_elements(1)
  m = LillyMol.MolFromSmiles("[Xx]C[Yy]")
  remove_all_non_natural_elements(m) == 2 || return false
  smiles(m) == "C" || return false
  LillyMol.set_auto_create_new_elements(0)
  true
end

function test_remove_explicit_hydrogens()::Bool
  m = LillyMol.MolFromSmiles("[H]C([H])([H])CC")
  valence_ok(m) || return false
  remove_explicit_hydrogens(m) == 3 || return false
  true
end

# Chop is not exported from base, not sure why it did not work.
# Seldom used...
function test_chop()::Bool
  m = LillyMol.MolFromSmiles("CN")
  LillyMol.chop(m, 1)
  smiles(m) == "C" || return false
  build_from_smiles(m, "CNN") || return false
  LillyMol.chop(m, 2)
  smiles(m) == "C" || return false
  true
end

function test_remove_bonds_to_atom()::Bool
  m = LillyMol.MolFromSmiles("FCN")
  remove_bonds_to_atom(m, 1)
  smiles(m) == "F.C.N" || return false
  true
end

function test_remove_bond()::Bool
  m = LillyMol.MolFromSmiles("OCN")
  remove_bond(m, 0)
  smiles(m) == "O.CN" || return false
  true
end

function test_remove_bond_between_atoms()::Bool
  m = LillyMol.MolFromSmiles("OCN")
  remove_bond_between_atoms(m, 0, 1)
  smiles(m) == "O.CN" || return false
  true
end

function test_remove_all_bonds()::Bool
  m = LillyMol.MolFromSmiles("CCC")
  remove_all_bonds(m)
  smiles(m) = "C.C.C" || return false
  true
end

function test_molecular_weight()::Bool
  m = LillyMol.MolFromSmiles("FC(Cl)(Br)CNC(=O)O")
  abs(220.4245 - molecular_weight(m)) < 0.001 || return false
  true
end

function test_molecular_weight_count_isotopes()::Bool
  m = LillyMol.MolFromSmiles("[1C]C")
  abs(molecular_weight_count_isotopes(m) - 16.03452) < 0.001 || return false
  true
end

function test_molecular_weight_ignore_isotopes()::Bool
  m = LillyMol.MolFromSmiles("[1C]C")
  abs(molecular_weight_ignore_isotopes(m) - 27.045221) < 0.001 || return false
  true
end

function test_highest_coordinate_dimensionality()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  highest_coordinate_dimensionality(m) == 0 || return false
  build_from_smiles(m, "C{{1,0,0}}")
  highest_coordinate_dimensionality(m) == 1 || return false
  build_from_smiles(m, "C{{1,1,0}}")
  highest_coordinate_dimensionality(m) == 2 || return false
  build_from_smiles(m, "C{{1,1,1}}")
  highest_coordinate_dimensionality(m) == 3 || return false
  true
end

function test_exact_mass()::Bool
  m = LillyMol.MolFromSmiles("CNOF")
  abs(exact_mass(m) - 65.027691992) < 0.001 || return false
  true
end

function test_number_fragments()::Bool
  m = Molecule()
  for i in 1:10
    smi = "C."^i
    smi *= "C"
    build_from_smiles(m, smi) || return false
    number_fragments(m) == (i + 1) || return false
  end
  true
end

function test_fragment_membership()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC.CCC")
  fragment_membership(m, 0) == 0 || return false
  fragment_membership(m, 1) == 1 || return false
  fragment_membership(m, 2) == 1 || return false
  fragment_membership(m, 3) == 2 || return false
  fragment_membership(m, 4) == 2 || return false
  fragment_membership(m, 5) == 2 || return false
  true
end

function test_fragment_membership_vector()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC.CCC")
  expected = [0, 1, 1, 2, 2, 2]
  fragm = Array{Int64}(undef, natoms(m))
  fragment_membership(m, fragm)
  for (i, f) in enumerate(fragm)
    expected[i] == f || return false
  end
  true
end

function test_atoms_in_fragment()::Bool
  smi = join(["C"^i for i in 1:10], '.')
# println("smile is $(smi)")
  m = Molecule()
  build_from_smiles(m, smi) || return false
  for i in 0:9
    atoms_in_fragment(m, i) == (i + 1) || return false
  end
  true
end

function test_get_atoms_in_fragment()::Bool
  m = Molecule()
  smi = join(["C"^i for i in 1:10], '.')
  build_from_smiles(m, smi)
  for i in 0:9
    s = get_atoms_in_fragment(m, i)
    length(s) == (i + 1) || return false
  end
  true
end

function test_largest_fragment()::Bool
  m = Molecule()
  smi = join(["C"^i for i in 1:10], '.')
  build_from_smiles(m, smi)
  largest_fragment(m) == 9 || return false
  true
end

function test_identify_spinach()::Bool
  m = LillyMol.MolFromSmiles("CC1C(O)C1N")
  spch = identify_spinach(m)
  spch[1] == 1 || return false
  spch[2] == 0 || return false
  spch[3] == 0 || return false
  spch[4] == 1 || return false
  spch[5] == 0 || return false
  spch[6] == 1 || return false
  true
end

function test_rings_in_fragment()::Bool
  m = LillyMol.MolFromSmiles("CCC.C1CC1.C1CC1C2CC2")
  rings_in_fragment(m, 0) == 0 || return false
  rings_in_fragment(m, 1) == 1 || return false
  rings_in_fragment(m, 2) == 2 || return false
  true
end

# Cannot figure out how to return a std::vector<Molecule>.
# This will change once I figure that out.
function test_create_components()::Bool
  m = LillyMol.MolFromSmiles("CCC.C1CC1.C1CC1C2CC2")
  nf = number_fragments(m)
  nf == 3 || return false
  frags = [Molecule() for i in 1:nf]
  create_components(m, frags)
  length(frags) == 3 || return false
  true
end

function test_create_subset()::Bool
  m = LillyMol.MolFromSmiles("Cc1ccccc1N")
  subset = [0, 1, 1, 1, 1, 1, 1, 0]
  length(subset) == natoms(m) || return false
  s = create_subset(m, subset)
  unique_smiles(s) == "c1ccccc1" || return false
  true
end

function test_create_subset_set_of_atoms()::Bool
  m = LillyMol.MolFromSmiles("Cc1ccccc1N")
  # Initializer list does not work.
  # subset = SetOfAtoms(1,2,3,4,5,6)
  subset = SetOfAtoms()
  for i in 1:6
    LillyMol.push!(subset, i)
  end
  s = create_subset(m, subset)
  unique_smiles(s) == "c1ccccc1" || return false
end

function test_reduce_to_largest_fragment()::Bool
  m = LillyMol.MolFromSmiles("CC")
  reduce_to_largest_fragment(m)
  smiles(m) == "CC" || return false
  build_from_smiles(m, "C.C")
  reduce_to_largest_fragment(m)
  smiles(m) == "C" || return false
  build_from_smiles(m, "C.N")
  reduce_to_largest_fragment(m)
  smiles(m) == "C" || return false
  build_from_smiles(m, "N.C")
  reduce_to_largest_fragment(m)
  smiles(m) == "N" || return false
  build_from_smiles(m, "C.CCCC")
  reduce_to_largest_fragment(m)
  smiles(m) == "CCCC" || return false
  build_from_smiles(m, "CCCC.C")
  reduce_to_largest_fragment(m)
  smiles(m) == "CCCC" || return false
  smi = join(Random.shuffle(["C"^i for i in 1:10]), '.')
  build_from_smiles(m, smi)
  reduce_to_largest_fragment(m)
  smiles(m) == "C"^10 || return false
  true
end

function test_reduce_to_largest_organic_fragment()::Bool
  m = LillyMol.MolFromSmiles("[Na].C")
  reduce_to_largest_organic_fragment(m)
  smiles(m) == "C" || return false
  true
end

function test_reduce_to_largest_fragment_carefully()::Bool
  m = LillyMol.MolFromSmiles("O=N(=O)c1cc(N(=O)=O)cc(N(=O)=O)c1O."* "OCN"^5)
  reduce_to_largest_fragment_carefully(m)
  smiles(m) == "OCN"^5 || return false
  true
end

function test_organic_only()::Bool
  m = LillyMol.MolFromSmiles("IC(Cl)(Br)C(F)NCOCSP")
  organic_only(m) || return false
  build_from_smiles(m, "[B]")
  organic_only(m) && return false
  build_from_smiles(m, "[Si]")
  organic_only(m) && return false
  build_from_smiles(m, "[Se]")
  organic_only(m) && return false
  true
end

function test_contains_non_periodic_table_elements()::Bool
  LillyMol.set_auto_create_new_elements(1)
  m = Molecule()
  build_from_smiles(m, "[He][Ll]O[W][Rl][D]") || return false
  contains_non_periodic_table_elements(m) || return false
  LillyMol.set_auto_create_new_elements(0)
  true
end

function test_longest_path()::Bool
  smi = "C"^64
  m = Molecule()
  build_from_smiles(m, smi) || return false
  longest_path(m) == (length(smi) - 1) || return false
  true
end

function test_atoms_between()::Bool
  m = Molecule()
  build_from_smiles(m, "CC1CC(C)C1") || return false
  btw = atoms_between(m, 0, 4)
  length(btw) == 3 || return false
  btw[0] == 1 || return false
  btw[1] == 2 || return false
  btw[2] == 3 || return false
  true
end

function test_implicit_hydrogens()::Bool
  m = LillyMol.MolFromSmiles("CCC(C)C(C)(C)C")
  implicit_hydrogens(m, 0) == 3 || return false
  implicit_hydrogens(m, 1) == 2 || return false
  implicit_hydrogens(m, 2) == 1 || return false
  implicit_hydrogens(m, 4) == 0 || return false
  build_from_smiles(m, "CC(=O)O") || return false
  implicit_hydrogens(m, 0) == 3 || return false
  implicit_hydrogens(m, 1) == 0 || return false
  implicit_hydrogens(m, 2) == 0 || return false
  implicit_hydrogens(m, 3) == 1 || return false
  build_from_smiles(m, "NCNCN(C)C#N") || return false
  implicit_hydrogens(m, 0) == 2 || return false
  implicit_hydrogens(m, 2) == 1 || return false
  implicit_hydrogens(m, 4) == 0 || return false
  implicit_hydrogens(m, 7) == 0 || return false
  build_from_smiles(m, "C[N+](C)C") || return false
  valence_ok(m) && return false
  implicit_hydrogens(m, 1) == 0 || return false
  build_from_smiles(m, "C[NH+](C)C") || return false
  valence_ok(m) || return false
  implicit_hydrogens(m, 1) == 1 || return false
  build_from_smiles(m, "SC") || return false
  implicit_hydrogens(m, 0) == 1 || return false
  true
end

function test_explicit_hydrogens()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  explicit_hydrogens(m, 0) == 0 || return false
  build_from_smiles(m, "[H]C") || return false
  explicit_hydrogens(m, 1) == 1 || return false
  build_from_smiles(m, "[H]C[H]") || return false
  explicit_hydrogens(m, 1) == 2 || return false
  build_from_smiles(m, "[H]C([H])[H]") || return false
  explicit_hydrogens(m, 1) == 3 || return false
  build_from_smiles(m, "[H]C([H])([H])[H]") || return false
  explicit_hydrogens(m, 1) == 4 || return false
  true
end

function test_hcount()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  hcount(m, 0) == 4 || return false
  build_from_smiles(m, "[H]C")
  hcount(m, 1) == 4 || return false
  build_from_smiles(m, "[H]C[H]")
  hcount(m, 1) == 4 || return false
  build_from_smiles(m, "[H]C([H])[H]")
  hcount(m, 1) == 4 || return false
  build_from_smiles(m, "[H]C([H])([H])[H]")
  hcount(m, 1) == 4 || return false
  true
end

function test_make_implicit_hydrogens_explicit()::Bool
  m = Molecule()
  build_from_smiles(m, "CC")
  make_implicit_hydrogens_explicit(m)
  smiles(m) == "C([H])([H])([H])C([H])([H])[H]" || return false
  true
end

function test_move_hydrogens_to_end_of_connection_table()::Bool
  m = Molecule()
  build_from_smiles(m, "CN")
  atomic_number(m, 0) == 6 || return false
  atomic_number(m, 1) == 7 || return false
  move_hydrogens_to_end_of_connection_table(m, 6)
  atomic_number(m, 0) == 7 || return false
  atomic_number(m, 1) == 6 || return false
  smiles(m) == "NC" || return false
  true
end

function test_pi_electrons()::Bool
  m = Molecule()
  build_from_smiles(m, "OCNC=NC#N")
  pi_electrons(m, 0) == 2 || return false
  pi_electrons(m, 1) == 0 || return false
  pi_electrons(m, 2) == 2 || return false
  pi_electrons(m, 3) == 1 || return false
  pi_electrons(m, 4) == 1 || return false
  pi_electrons(m, 5) == 2 || return false
  pi_electrons(m, 6) == 2 || return false
  true
end

function test_lone_pair_count()::Bool
  m = Molecule()
  build_from_smiles(m, "OCNC=NC#N")
  lone_pair_count(m, 0) == 2 || return false
  lone_pair_count(m, 1) == 0 || return false
  lone_pair_count(m, 2) == 1 || return false
  lone_pair_count(m, 3) == 0 || return false
  lone_pair_count(m, 4) == 1 || return false
  lone_pair_count(m, 5) == 0 || return false
  lone_pair_count(m, 6) == 1 || return false
  true
end

function test_aromatic_atom_count()::Bool
  m = LillyMol.MolFromSmiles("Cc1nccc1")
  aromatic_atom_count(m) == 5 || return false
  true
end

function test_aromatic_ring_count()::Bool
  m = LillyMol.MolFromSmiles("c12ccccc2ccc(c1)c1ccccc1")
  aromatic_ring_count(m) == 3 || return false
  true
end

function test_saturated()::Bool
  m = LillyMol.MolFromSmiles("CC=Cc1ccccc1")
  saturated(m, 0) || return false
  saturated(m, 1) && return false
  saturated(m, 2) && return false
  saturated(m, 3) && return false
  true
end

function test_chiral_centres()::Bool
  m = Molecule()
  build_from_smiles(m, "C[C@H](F)N") || return false
  chiral_centres(m) == 1 || return false
  true
end

function test_chiral_centre_at_atom()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N")
  c = chiral_centre_at_atom(m, 1)
  top_front(c) == 0 || return false
  left_down(c) == 2 || return false
  right_down(c) == 3 || return false
  # TODO:ianwatson figure out lone pairs and implicit hydrogens
  true
end

function test_chiral_centre_in_molecule_not_indexed_by_atom_number()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N[C@@H](C)CC")
  c = chiral_centre_in_molecule_not_indexed_by_atom_number(m, 0)
  centre(c) == 1 || return false
  c = chiral_centre_in_molecule_not_indexed_by_atom_number(m, 1)
  centre(c) == 4 || return false
  true
end

function test_remove_chiral_centre_at_atom()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N[C@@H](C)CC")
  remove_chiral_centre_at_atom(m, 4)
  smiles(m) == "C[C@H](F)NC(C)CC" || return false
  true
end

function test_remove_all_chiral_centres()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N[C@@H](C)CC")
  remove_all_chiral_centres(m)
  smiles(m) == "CC(F)NC(C)CC" || return false
  true
end

function test_invert_chirality_on_atom()::Bool
  m = LillyMol.MolFromSmiles("C[C@H](F)N[C@@H](C)CC")
  invert_chirality_on_atom(m, 4)
  smiles(m) == "C[C@H](F)N[C@H](C)CC" || return false
  true
end

function test_smarts_equivalent_for_atom()::Bool
  m = LillyMol.MolFromSmiles("Cc1ncccc1O")
  smarts_equivalent_for_atom(m, 0) == "[CD1H3v1]" || return false
  smarts_equivalent_for_atom(m, 1) == "[cD3H0v4;r6]" || return false
  true
end

function test_smarts()::Bool
  m = LillyMol.MolFromSmiles("Cn1cnc2n(C)c(=O)n(C)c(=O)c12")
  smarts(m) == "[C][n]1[c][n][c]2[c]1[c](=[O])[n]([C])[c](=[O])[n]2[C]" || return false
  true
end

function test_atom_map_number()::Bool
  m = Molecule()
  build_from_smiles(m, "[CH3:1]C[CH3:3]")
  atom_map_number(m, 0) == 1 || return false
  atom_map_number(m, 1) == 0 || return false
  atom_map_number(m, 2) == 3 || return false
  true
end

function test_set_atom_map_number()::Bool
  m = LillyMol.MolFromSmiles("CCC")
  set_atom_map_number(m, 1, 8)
  smiles(m) == "C[CH2:8]C" || return false
  atom_map_number(m, 1) == 8 || return false
  set_atom_map_number(m, 1, 0)
  smiles(m) == "CCC" || return false
  true
end

function test_set_include_atom_map_with_smiles()::Bool
  m = Molecule();
  build_from_smiles(m, "[CH3:7]C") || return false
  set_include_atom_map_with_smiles(false)
  atom_map_number(m, 0) == 7 || return false
  # This is kind of unfortunate, perhaps this could be fixed/changed.
  smiles(m) == "[CH3]C" || return false
  set_include_atom_map_with_smiles(true)
  true
end

function test_atom_with_atom_map_number()::Bool
  m = Molecule()
  build_from_smiles(m, "C[CH2:4]C")
  atom_with_atom_map_number(m, 4) == 1 || return false
  atom_with_atom_map_number(m, 7) < 0 || return false
  true
end

function test_reset_all_atom_map_numbers()::Bool
  m = LillyMol.MolFromSmiles("C[CH2:5]C")
  reset_all_atom_map_numbers(m)
  smiles(m) == "CCC" || return false
  true
end

function test_unset_unnecessary_implicit_hydrogens_known_values()::Bool
  m = Molecule()
  build_from_smiles(m, "[OH][CH2][NH2]")
  unset_unnecessary_implicit_hydrogens_known_values(m)
  smiles(m) == "OCN" || return false
  true
end

function test_discern_chirality_from_3d_structure()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,1,1}}C{{0,0,0}}(N{{0,1,-1}})(O{{-1,-1,0}})C{{1,-1,0}}C") || return false
  smiles(m) == "CC(N)(O)CC" || return false
  discern_chirality_from_3d_structure(m)
  smiles(m) == "C[C@](N)(O)CC" || return false
  setz(m, 0, -1.0)
  setz(m, 2, 1.0)
  discern_chirality_from_3d_structure(m)
  smiles(m) == "C[C@@](N)(O)CC" || return false
  true
end

function test_bonds_between()::Bool
  m = LillyMol.MolFromSmiles("Oc1ccc(C)cc1")
  bonds_between(m, 0, 1) == 1 || return false
  bonds_between(m, 0, 2) == 2 || return false
  bonds_between(m, 0, 3) == 3 || return false
  bonds_between(m, 0, 4) == 4 || return false
  bonds_between(m, 0, 5) == 5 || return false
  bonds_between(m, 0, 6) == 3 || return false
  bonds_between(m, 0, 7) == 2 || return false
  bonds_between(m, 1, 2) == 1 || return false
  bonds_between(m, 1, 3) == 2 || return false
  bonds_between(m, 1, 4) == 3 || return false
  bonds_between(m, 1, 5) == 4 || return false
  bonds_between(m, 1, 6) == 2 || return false
  bonds_between(m, 1, 7) == 1 || return false
  true
end

function test_returns_vector()::Bool
  v = LillyMol.returns_vector(1, 2)
  v == [1, 2] || return false
  true
end

function test_standardise()
  standardise = ChemicalStandardisation()
  activate_all(standardise)
  m = Molecule()
  build_from_smiles(m, "CC[N+](=O)[O-]")
  process(standardise, m)
  return smiles(m) == "CCN(=O)=O"
end

#for (index,atom) in enumerate(m)
#  print("atom $(index) type $(atomic_number(atom)) connected to")
#  for bond in atom
#    nbr = other(bond, index)
#    print(" $(nbr)")
#  end
#  println("")
#end

boobar()
@test test_empty_molecule()
@test test_set_of_atoms()
@test test_build_from_smiles()
@test test_set_name()
@test test_natoms()
@test test_natoms_atomic_number()
@test test_natoms_element()
@test test_length_molecule()
@test test_atomic_number()
@test test_atomic_number_in()
@test test_nedges()
@test test_molecular_formula()
@test test_valence_ok()
@test test_standardise()
@test test_nrings()
@test test_nrings_atom()
@test test_nrings_size()
@test test_is_ring_atom()
@test test_fused_system_size()
@test test_fused_system_identifier()
@test test_rings_with_fused_system_identifier()
@test test_in_same_ring()
@test test_in_same_aromatic_ring()
@test test_in_same_ring_system()
@test test_ring_membership()
@test test_rings_containing_both()
@test test_is_part_of_fused_ring_system()
@test test_ring()
@test test_rings()
@test test_ring_containing_atom()
@test test_label_atoms_by_ring_system()
@test test_label_atoms_by_ring_system_including_spiro_fused()
@test test_nrings_including_non_sssr_rings()
@test test_non_sssr_rings()
@test test_non_sssr_ring()
@test test_is_spiro_fused()
@test test_is_halogen()
@test test_ncon_molecule()
@test test_nbonds_molecule()
@test test_maximum_connectivity()
@test test_connections_molecule()
@test test_isotopically_labelled_smiles()
@test test_is_aromatic()
@test test_getindex_molecule()
@test test_number_formally_charged_atoms()
@test test_net_formal_charge()
@test test_bond_molecule()
@test test_bond_between_atoms()
@test test_number_symmetry_classes()
@test test_symmetry_class()
@test test_symmetry_equivalents() broken=true
@test test_symmetry_classes()
@test test_attached_heteroatom_count()
@test test_bond_length()
@test test_bond_angle()
#@test test_dihedral_angle()
#@test test_signed_dihedral_angle()
#@test test_set_bond_length()
#@test test_set_bond_angle()
#@test test_set_dihedral_angle()
#@test test_xyx_molecule()
#@test test_distance_between_atoms()
#@test tests_longest_intra_molecular_distance()
@test test_add_bond()
@test test_are_bonded()
@test test_bond_between_atoms()
@test test_bond_list()
@test test_add_molecule()
@test test_remove_atom()
@test test_remove_atoms()
@test test_remove_atoms_vector()
@test test_delete_fragment()
@test test_remove_fragment_containing_atom()
@test test_remove_all()
@test test_set_auto_create_new_elements()
@test test_set_atomic_symbols_can_have_arbitrary_length()
@test test_remove_all_non_natural_elements()
@test test_remove_explicit_hydrogens()
@test test_chop()
@test test_remove_bonds_to_atom()
@test test_remove_bond()
@test test_remove_bond_between_atoms()
@test test_remove_all_bonds()
@test test_molecular_weight()
@test test_molecular_weight_count_isotopes()
@test test_molecular_weight_ignore_isotopes()
@test test_highest_coordinate_dimensionality()
@test test_exact_mass()
@test test_number_fragments()
@test test_fragment_membership()
@test test_fragment_membership_vector()
@test test_atoms_in_fragment()
@test test_get_atoms_in_fragment()
@test test_largest_fragment()
@test test_identify_spinach()
@test test_rings_in_fragment()
@test test_create_components() skip=true
@test test_returns_vector()
@test test_create_subset()
@test test_create_subset_set_of_atoms()
@test test_reduce_to_largest_fragment()
@test test_reduce_to_largest_organic_fragment()
@test test_reduce_to_largest_fragment_carefully()
@test test_organic_only()
@test test_contains_non_periodic_table_elements()
@test test_longest_path()
@test test_bonds_between()
@test test_atoms_between()
@test test_implicit_hydrogens()
@test test_explicit_hydrogens()
@test test_hcount()
@test test_make_implicit_hydrogens_explicit()
@test test_move_hydrogens_to_end_of_connection_table()
@test test_pi_electrons()
@test test_lone_pair_count()
@test test_aromatic_atom_count()
@test test_aromatic_ring_count()
@test test_saturated()
@test test_chiral_centres()
@test test_chiral_centre_at_atom()
@test test_chiral_centre_in_molecule_not_indexed_by_atom_number()
@test test_remove_chiral_centre_at_atom()
@test test_remove_all_chiral_centres()
@test test_invert_chirality_on_atom()
@test test_smarts_equivalent_for_atom()
@test test_smarts()
@test test_atom_map_number()
@test test_set_atom_map_number()
@test test_set_include_atom_map_with_smiles()
@test test_atom_with_atom_map_number()
@test test_reset_all_atom_map_numbers()
@test test_unset_unnecessary_implicit_hydrogens_known_values()
@test test_discern_chirality_from_3d_structure()
