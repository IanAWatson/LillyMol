# LillyMol Python
This release contains experimental python bindings for some parts of LillyMol.

We find that using LillyMol from python can be a very effective means of prototyping
an idea, or doing things that are not time sensitive.

This first release comprises three main components

1. The Molecule object
2. Substructure Searching
3. Reaction Handling

More functionality will become available.

The current python bindings may not reflect final names or functionality, this
is a work in progress, but has already proven useful.

## Background
Python bindings for LillyMol were implemented at Google in 2019 using 
[Clif](https://google.github.io/clif/clif/python/). While this appeared to work
well, and Clif was easy to use inside Google, it proved to be very difficult
outside Google. Instead [pybind11](https://github.com/pybind/pybind11) has been
used here. It is an amazing template metaprogramming tour-de-force, which
ultimately seems to work well.

Unfortunately there do appear to be many instances of needless copies happening
between C++ and Python. Perhaps these could be lessened via careful inspection,
but for now, there is no claim that this is as fast as things could be.

## Building and installation
Your python environment *must* include pybind11. Normally
```
pip install pybind11
```
will accomplish this.

Normally the python bindings are built as part of the default build,
the script [build_linux.sh](/src/build_linux.sh), but if
you wish to compile separately that can be done via
```
bazelisk --output_user_root=/local/disk/ian build --cxxopt=-DTODAY="$(date +%Y-%b-%d)" --cxxopt=-DGIT_HASH="$(git rev-parse --short --verify HEAD)" --local_cpu_resources=10 -c opt pybind:all
```
This generates several `*.so` files in `bazel-bin/pybind`. In addition,
LillyMol now has several run-time dependencies, and these also need to be made
available.

### Local development with `run_python.sh`

For local development, the script `copy_shared_libraries.sh` in the `src`
directory copies the pybind extension modules and their LillyMol shared-library
dependencies out of `bazel-bin` and into `${LILLYMOL_HOME}/lib`.

```
cd ${LILLYMOL_HOME}/src
bazel build //pybind:all
./copy_shared_libraries.sh ../lib
```

Once the shared libraries are copied to `${LILLYMOL_HOME}/lib`, the top-level
`run_python.sh` script can be used to invoke Python with `PYTHONPATH` and
`LD_LIBRARY_PATH` set appropriately.

```
${LILLYMOL_HOME}/run_python.sh my_script.py
```

This is the simplest workflow while developing LillyMol itself.

### Building a wheel

A wheel can be built from the same Bazel-generated pybind modules. This does not
compile C++ via setuptools; Bazel still does the C++ build. The wheel machinery
just packages the prebuilt extension modules and their private LillyMol shared
libraries.

```
cd ${LILLYMOL_HOME}/src
bazel build //pybind:all
./copy_shared_libraries.sh ../lib

cd ${LILLYMOL_HOME}/python
./scripts/stage_wheel_files.sh
python -m pip install build wheel setuptools
python -m build --wheel
```

The wheel is written to `${LILLYMOL_HOME}/python/dist`. The staging script
creates `${LILLYMOL_HOME}/python/prebuilt`, which is generated data and should
not be committed.

The wheel includes the pybind extension modules, selected generated protobuf
Python modules, and private LillyMol shared libraries. On Linux the staged shared
objects are patched with an `$ORIGIN` runpath so installed extension modules can
find the private shared libraries shipped beside them, without requiring
`LD_LIBRARY_PATH`.

The wheel can be smoke-tested in a clean target directory with something like

```
rm -rf /tmp/lillymol_wheel_test
python -m pip install --no-deps --target /tmp/lillymol_wheel_test \
  ${LILLYMOL_HOME}/python/dist/lillymol-*.whl

env -u LD_LIBRARY_PATH PYTHONPATH=/tmp/lillymol_wheel_test python - <<'PYTHON_SMOKE_TEST'
import lillymol
import lillymol_query
import lillymol_tools
m = lillymol.MolFromSmiles("CCO ethanol")
q = lillymol_query.QueryFromSmarts("[OD1,OD2]")
assert q in m
ctx = lillymol_tools.GFPContext.standard()
fp = ctx.fingerprint(m)
assert abs(ctx.distance(fp, fp)) < 1.0e-6
print("LillyMol wheel smoke test ok")
PYTHON_SMOKE_TEST
```

See `${LILLYMOL_HOME}/python/README.md` for the wheel staging details.

See `MODULE.bazel` for how we configured the local python and pybind11 installs.

## Philosophy
LillyMol has no concept of changeable and unchangeable molecules. Any molecule
can be altered at any time. This simplifies interaction with Molecules. In the
python implementation it does raise some risks of errors, while making certain operations
easier. Read on...

This works because the LillyMol Molecule is a very lazy object. It never computes
things like fragment membership, ring membership, aromaticity or canononical
ordering unless requested. So if you remove an atom, or bond, it will destroy any
information it has about those derived quantities. Only if requested will any
be recomputed.

This has many advantages. For example if a molecule is built from a smiles
and then the only thing ever requested is the number of atoms, that will be
very cheap. If the number of fragments is requested, then only fragment membership
is computed. Even if the most expensive derived property, the canonical order,
is requested, the actual unique smiles will not be generated unless requested.

While this generally works well, there is one caveat. Because of this, things
like Rings and Bonds, by default, do not know if they are aromatic or not, or
if they are in a ring or not. Since neither one knows anything about being in
a Molecule, the following will fail
```
benzene.build_from_smiles('c1ccccc1')
benzene.ring(0).is_aromatic()    # Is the first ring aromatic?
```
But this will work
```
if benzene.is_aromatic(0):    # Is the first atom aromatic?
  ....
benzene.ring(0).is_aromatic()  # Is the first ring aromatic?
```
In this case `benzene.is_aromatic(0)` meant that the Molecule needed to 
compute fragment membership, ring membership and aromaticity. Then
when the first ring, `ring(0)` was queried, it now knew that it was
aromatic.

Again, this is only an issue if you will be querying molecular properties
via Atoms, Bonds or Rings. If you get this information by asking the
Molecule, this happens automatically. If you are going to query individual
objects for things like aromaticity and ring membership start with
```
mol.compute_aromaticity_if_needed()
```
which will ensure that aromaticity information is propagated throughout
all parts of the molecule - and if that has already been done, this does
nothing.

Atoms do not know if they are in a ring or not, nor whether 
they are aromatic, or what fragment they are in. Fragment membership, 
ring  membership and aromaticity are molecular properties.
Same with chirality. Bonds are slightly different, in that Bonds do
know their ring membership, and if they are aromatic or not. But again
these quantities are not computed by default.

## Using LillyMol Python
The Molecule functionality must be imported
```
from lillymol import *
```
For those familiar with RDKit, this enables
```
m = MolFromSmiles('c1ccccc1')
```
and if an invalid smiles is encountered, None will be returned. There
is also a list form of this
```
mols = MolFromSmiles(["C", "CC", "C1CC1"])
```
which returns a list of molecules. This may offer speed advantages
depending on the structure of the program.

And for clarity
```
mol = LillyMolFromSmiles("C methane")
```
also works.


There are other means by which molecules can enter the system.

```
m = Molecule()
```
instantiates an empty molecule. That can then be constructed via something like
```
m.build_from_smiles('c1ccccc1')
```
which returns True if successful. The same Molecule can be re-used any number of
times with the `build_from_smiles` method, which first cleans out all atoms before
starting. See discussion of the addition operators later.

Molecules can also be read from files, which will likely be the most common case.

```
from lillymol_io import *

reader = Reader()
if not reader.open('/path/to/file.smi', FileType.SMI):
  logging.error('Cannot open...')

for m in reader:
  print(m)
```
This will print the smiles and name of each molecule. If reading a 
.sdf file, use `FileType.SDF`. LillyMol has a wide variety of directives
for reading .sdf files, those need to be made available via Python.

Note that there will never be a None molecule returned. If a connection
table error is encountered, reading will cease. The Reader class has
a 'set_set_connection_table_errors_allowed' method, which allows you
to set the number of otherwise fatal errors that are ignored. Warnings
will flash by on stderr, but nothing will show up in python. 

Lillymol has always operated on the principle that your input should
be correct. That said, it would not be hard to add an option to return
a None molecule in the event of an otherwise ignored error.

Note that in LillyMol a molecule with an invalid valence is a valid molecule.
If you don't want molecules with valence errors, use the valence_ok() method
to check each molecule and skip those having an invalid valence.

A more pythonic way of reading structures is available as
```
  with ReaderContext('/path/to/file.sdf') as reader:
    for mol in reader:
      for atom in mol:
        ...
```
This last example shows that a molecule is iterable, and a stream
of Atom objects is returned. In LillyMol, Atoms know nothing about
a Molecule, so if you want the atom number something like
```
for atom_number, atom in enumerate(mol):
  print(f'atom {atom_number} atomic number {atom.atomic_number()}')
```
the atoms are iterated in sequential order, so enumerate works.

While an Atom does not know anything about being part of a Molecule,
the Molecule can query each of its atoms for properties associated
with each atom. Therefore
```
m.atomic_number(3)
m[3].atomic_number()
```
both report the atomic number of atom 3. The first queried the
molecule, and the second retrieved the third atom and asked it
for the atomic number. In terms of efficiency, the first method
will be more efficient in python.

And for the greatest simplicity in getting a list of molecules
into python
```
mols = slurp(fname)
```
will read all molecules from `fname` into a list. Note that this
only works if `fname` has the proper suffix, and it will not work
if trying to read stdin. Returns None if anything goes wrong. Note
that if any molecule fails nothing is returned.

## Molecule Methods
The most common methods for a Molecule currently implemented are

| Method | Description |
| ------ | ----------- |
| name() | The name of the molecule |
| set_name(string) | Set the name |
| natoms() | Number of atoms (explicit atoms only) |
| empty() | True if there are no atoms in the molecule |
| GetNumAtoms() | Number of atoms (explicit atoms only) |
| resize(n) | keep only the first n atoms |
| nedges() | Number of bonds |
| bonds() | Iterable collection of Bonds |
| nrings() | Number of SSSR rings |
| nrings(atom) | Ring membership of 'atom' |
| is_ring_atom(atom) | True if 'atom' is in a ring |
| IsInRing(atom) | True if 'atom' is in a ring |
| in_ring_of_given_size(atom, rsize) | True if 'atom' is in a ring of size 'rsize' |
| IsAtomInRingOfSize(atom, rsize) | True if 'atom' is in a ring of size 'rsize' |
| ring_bond_count(atom) | Number of ring bonds involving 'atom' |
| get_ring_membership() | Ring membership for each atom |
| number_ring_systems() | Number of ring systems in the molecule |
| fused_system_identifier(atom) | Fused system identifier for 'atom' |
| fused_system_size(atom) | Size of fused system containing 'atom' |
| label_atoms_by_ring_system() | Fused system identifier for each atom |
| label_atoms_by_ring_system_including_spiro_fused() | Ring systems span spiro fusions |
| ring(number) | Fetch a particular ring |
| rings() | Iterable collection of all rings |
| in_same_ring(a1, a2) | True if a1 and a2 are in the same ring | 
| in_same_ring_system(a1, a2) | True if a1 and a2 are in the same ring system |
| largest_ring_size() | Number of atoms in largest ring |
| number_ring_systems() | Number of ring systems - naphthalene counts as 1 |
| is_spiro_fused(atom) | True if 'atom' is a spiro fusion |
| amw() | Average Molecular Weight |
| molecular_formula() | Molecular Formula |
| natoms(atomic_number) | Number of atoms with atomic number |
| natoms(atomic_symbol) | Number of atoms with atomic symbol |
| exact_mass() | Exact Mass |
| ncon(atom) | Number of edges (bonds) to 'atom' |
| connections(atom) | List of all atoms attached to 'atom' |
| other_atom(atom, n) | Atom number of n'th connection to 'atom' |
| attached_heteroatom_count(atom) | Number of heteroatoms attached to 'atom ' |
| is_aromatic(atom) | True if 'atom' is aromatic |
| compute_aromaticity_if_needed() | Compute fragments, rings and aromaticity |
| aromatic_atom_count() | Number of aromatic atoms |
| aromatic_ring_count() | Number of aromatic rings |
| atomic_number(atom) | Atomic number of 'atom' |
| atomic_symbol(atom) | Atomic symbol of 'atom' |
| set_atomic_number(atom, atomic_number) | Change an element |
| add_bond(atom1, atom2, btype) | Add a bond |
| set_bond_type_between_atoms(atom1, atom2, btype) | Change an existing bond |
| is_halogen(atom) | True if 'atom' is a Halogen |
| remove_atom(atom) | Remove an atom |
| remove_atoms(list, flag) | Remove all atoms where list[i] == flag |
| remove_atoms(Set_of_Atoms) | Remove the atoms in the set |
| remove_atoms(numpy_array, flag) | Remove atoms where numpy_array[i] == flag |
| remove_non_periodic_table_elements() | Remove any non-natural atoms |
| remove_all(atomic_number) | Remove all atoms with atomic_number |
| move_to_end_of_connection_table(z) | Move all atoms with atomic number z to end of connection table |
| chop(n) | Remove the last 'n' atoms in the molecule |
| organic_only() | True if only C, N, O, F, P, S, Cl, Br, I |
| remove_explicit_hydrogens() | Remove explicit Hydrogens |
| RemoveHs() | Remove explicit Hydrogens |
| implicit_hydrogens(atom) | Number of implicit Hydrogens on 'atom' |
| explicit_hydrogens(atom) | Number of explicit Hydrogens attached to 'atom' |
| make_implicit_hydrogens_explicit() | Make implicit Hydrogens into explicit Atoms |
| AddHs() | Make implicit Hydrogens into explicit Atoms |
| hcount(atom) | Sum of implicit and explicit Hydrogens for 'atom' |
| implicit_hydrogens_known(atom) | True if 'atom' has [] in smiles |
| saturated(atom) | True if 'atom' is fully saturated |
| pi_electrons(atom) | Pi electrons on 'atom' |
| lone_pair_count(atom) | Lone pairs on 'atom' |
| remove_all(atomic_number) | Remove all instances of that atom type |
| remove_bonds_to_atom(atom) | Remove all bonds involving 'atom' |
| remove_edge(edge) | Remove a bond by bond number |
| remove_bond_between_atoms(a1, a2) | Remove bond between atoms |
| remove_all_bonds() | All atoms become their own fragment |
| smarts_equivalent_for_atom(atom) | Smarts for 'atom' |
| number_fragments() | number of fragments |
| fragment_membership(atom) | Fragment number for 'atom' |
| atoms_in_fragment(frag) | Number of atoms in a fragment |
| delete_fragment(frag) | Delete a fragment |
| remove_fragment_containing_atom(atom) | Remove fragment containing atom |
| reduce_to_largest_fragment() | Discard all but the largest fragment |
| reduce_to_largest_fragment_carefully() | Contains heuristics to do a better selection |
| get_fragment_membership() | Return a list of fragment memberships |
| create_components() | Return a list of Molecules from a multi fragment molecule |
| to_scaffold() | Remove all non-scaffold atoms |
| canonical_rank(atom) | Canonical rank for 'atom' |
| canonical_ranks() | Canonical rank of each atom |
| symmetry_class(atom) | Symmetry class for 'atom' |
| symmetry_equivalents(atom) | A list of the atoms equivalent to 'atom' |
| number_symmetry_classes() | Number symmetry classes |
| build_from_smiles(smi) | Build from smiles |
| smiles() | Smiles |
| unique_smiles() | Unique smiles |
| random_smiles() | Random smiles |
| isotopically_labelled_smiles() | Each atom has isotope according to atom number |
| unique_kekule_smiles() | Unique Kekule form (expensive to compute) |
| smarts() | Molecule as smarts - does not work for searching |
| smiles_starting_with_atom(atom) | smiles where 'atom' is the first atom |
| smiles_atom_order() | atom order in must recent smiles produced |
| are_bonded(a1, a2) | True if a1 and a2 are bonded |
| add(Molecule other) | Add the atoms and bonds of 'other' |
| remove_hydrogens_known_flag_to_fix_valence_errors | Remove problematic square brackets |
| unset_unnecessary_implicit_hydrogens_known_values() | Try to fix certain valence problems |
| formal_charge(atom) | Formal charge on atom |
| set_formal_charge(atom) | Set formal charge on atom |
| has_formal_charges() | True if any atom has a formal charge |
| number_formal_charges() | Number of formally charged atoms |
| net_formal_charge() | Net formal charge |
| number_chiral_centres() | Number of chiral centres |
| remove_all_chiral_centres() | Remove all chiral centres |
| chiral_centre(atom) | Return the Chiral_Centre on 'atom' |
| invert_chirality_on_atom(atom) | Invert chirality |
| chiral_centres() | Iterable list of Chiral_Centre |
| isotope(atom) | Isotope on 'atom' |
| set_isotope(atom, iso) | Set isotope |
| set_isotopes(Set_of_Atoms, iso) | Set isotope for atoms in the set |
| set_isotopes(numpy_array) | Set each isotope |
| remove_isotopes() | Remove all isotopes |
| number_isotopic_atoms() | Number of atoms with non zero isotopes |
| bonds_between(a1, a2) | Bonds between atoms |
| longest_path() | Longest through bond path |
| atoms_on_shortest_path(a1, a2) | Set_of_Atoms holding atoms on shortest path between a1 and a2 |
| down_the_bond(a1, a2) | Return all atoms found by looking down the a1->a2 bond. May return None |
| atoms_by_radius(starting_atoms, max_radius) | List of Set_of_Atoms grouped by exact minimum bond distance |
| atom_map_number(atom) | Atom map number on 'atom' |
| set_atom_map_number(atom, nbr) | Set atom map number |
| reset_atom_map_numbers() | Remove all atom map numbers |
| atom_with_atom_map_number(number) | Atom with atom map number |
| bond_length(a1, a2) | Bond distance |
| bond_angle(a1, a2, a3) | Bond angle |
| dihedral_angle(a1, a2, a3, a4) | Dihedral angle |
| signed_dihedral_angle(a1, a2, a3, a4) | Dihedral angle, may be negative |
| distance_between_atoms(a1, a2) | Distance between atoms - bonded or not |
| longest_intra_molecular_distance() | Longest inter atom distance |
| x(atom) | x coordinate of 'atom' |
| y(atom) | y coordinate of 'atom' |
| z(atom) | z coordinate of 'atom' |
| setx(atom, x) | Set x coordinate of 'atom' |
| sety(atom, y) | Set y coordinate of 'atom' |
| setz(atom, z) | Set z coordinate of 'atom' |
| setxyz(atom, x, y, z) | Set coordinates of 'atom' |
| translate(x, y, z) | Translate atoms |
| highest_coordinate_dimensionality() | Will be 3 of 3D coordinates available |
| discern_chirality_from_3d_structure() | Use geometry to discern chiral centres |
| dihedral_scan(atom, atom, angle, bump_check | return list of coordinate sets |
| non_sssr_rings() | Number of non Smallest Set of Smallest Rings rings |
| non_sssr_ring(i) | The i'th non-SSSR ring |
| has_partial_charges() | True if the molecule has partial charges |
| invalidate_charges() | Discard any partial charge information stored |
| partial_charge_type() | The kind of partial charges stored |
| compute_Abraham_partial_charges() | Abraham partial charges |
| compute_Gasteiger_partial_charges() | Gasteiger partial charges |
| compute_Huckel_partial_charges() | Huckel partial charges |
| compute_Gasteiger_Huckel_partial_charges() | Gasteiger Huckel partial charges |
| compute_Del_Re_partial_charges() | Del Re partial charges |
| compute_Pullman_partial_charges() | Del Re partial charges |
| \__eq__ | True if m1 == m2. Will use unique smiles if necessary |
| m1 += m2 | Adds atoms and bonds from m2 to m1 |
| m1 + m2 | Returns a new molecule containing m1 and m2 |
| \__iter__ | List of Atoms |
| \__getitem__ | Get i'th atom |
| \__len__ | Number of atoms |
| \__eq__ | True if molecules contain same structures |
| \__contains__ | True of molecule contains atomic number |
| valence_ok() | True if all atoms have an OK valence |
| ok | True if the internal state of the Molecule is ok |
| debug_string() | String representation of internal state: print(m.debug_string())|
| ----- | ----- |

## Atom Methods
As mentioned previously, LillyMol Atoms are faily simple, and have no idea
that they are part of a Molecule. The only atributes an atom has is

* Pointer to an Element
* Isotope
* Coordinates
* Implicit Hydrogen info
* Formal Charge
* Number of bonds
* List of Bonds involving the Atom.
* Atom Map number

and some more obscure things.

The Atom object supports
| Method | Description |
| ------ | ----------- |
| atomic_number() | Atomic number |
| atomic_weight() | Atomic weight |
| exact_mass() | exact_mass |
| ncon() | Number of connections/edges |
| nbonds() | Number of bonds, single=1, double=2, triple=3 |
| formal_charge() | Formal charge |
| other(atom_number, ndx) | atom number of the 'ndx' connection |
| is_bonded_to(atom) | True if atom is bonded to 'atom' |
| valence_ok(atom) | True if valence ok |
| fully_saturated() | True if nbonds() == ncon() |
| atom_map() | atom map number |
| connections(atom) | iterable list of atoms attached |
| implicit_hydrogens | number of implicit hydrogens attached |
| \__iter__ | List of Bonds attached |
| \__contains__ | True if atom is bonded to |
| \__len__ | Number of connections |

In additon an Atom object inherits from an object that holds coordinates. Subsequent
versions will enable more of that functionality. For now the subtraction operator
returns the distance between two atoms, although long term this must be changed
so that subtraction of two atoms returns the vector between them.
```
m.build_from_smiles("C{{0,0,0}}C{{1,1,1}}"))
m[0] - m[1]
```
reports sqrt(3). For now...

A common construct might be (count the number of carbon=,#nitrogen bonds)
```
  result = 0
  for i,atom in enumerate(mol):
    if not atom.atomic_number() == 7:
      continue
    for bond in atom:
      if b.is_single_bond():
        continue
      other = b.other(i)
      if m.atomic_number(other) == 6:
        result += 1
```
But as is often the case, there is a more efficient way of doing this. The
above will visit each Bond twice - since each atom knows about all Bonds.
Traversing the bond list results in each Bond being examined only once.
```
  result = 0
  for bond in m.bond_list():
    if b.is_single_bond():
      continue
    a1 = b.a1()
    a2 = b.a2()
    if m.atomic_number(a1) == 6 && m.atomic_number(a2) == 7:
      result += 1
    elif m.atomic_number(a1) == 7 && m.atomic_number(a2) == 6:
      result += 1
```
Knowing when to solve a problem by traversing atoms and when to traverse
bonds can be hard.

## Bond Methods
Again, the Bond class really does not know much.

* The two atoms that define the bond.
* The bond type
* Ring membership
* Directional or not

And a couple of other things.

| Method | Description |
| ------ | ----------- |
| a1() | The first  atom |
| a2() | The second atom |
| btype() | The bond type |
| other(atom_number) | The atom that is not 'atom_number' |
| involves(atom_number) | True if 'atom_number' is a1 or a2 |
| is_single_bond() | True if the bond is a single bond |
| is_double_bond() | True if the bond is a double bond |
| is_triple_bond() | True if the bond is a triple bond |
| is_aromatic() | True if the bond is aromatic |
| nrings() | The number of rings involving this bond |
| is_directional(() | True if bond is directional |
| IsInRing() | True if bond is in a ring |
| GetBeginAtomIdx() | Same as a1() |
| GetEndAtomIdx() | Same as a2() |
| GetBondType() | Same as btype() |
| \__contains__ | involves() |

## Set_of_Atoms Methods
Set_of_Atoms objects are used extensively in LillyMol. Despite the name,
they are actually just vectors of atom numbers, with no requirement for
uniqueness. That said, it would be a very unusual application where a
Set_of_Atoms contained duplicate atom numbers.

| Method | Description |
| ------ | ----------- |
| empty() | True of the set is empty |
| size() | Number of items |
| scatter(list, value) | Set values to 'value' |
| \__len__ | Number of items |
| \__getitem__ | Access via [i] |
| \__iter__ | Access atoms via iterators |
| \__contains__ | Is atom included |

The C++ version contains several gather and scatter type methods. Other methods may
be added.

## Ring Methods
A Molecule may have rings. Ring's are just Set_of_Atoms's that have some
extra information.

* Aromaticity
* Fused Ring neighbours.

In addition, the atoms in the Ring are ordered, so if they are iterated, they will
trace out a bonded path through the ring.

| Method | Description |
| ------ | ----------- |
| size() | size |
| ring_number() | unique ring number |
| fragment_membership() | fragment number containing ring |
| fused_system_identifier() | fused system number containing this ring |
| is_fused() | True if ring is fused to another ring |
| fused_ring_neighbours() | Number of fused neighbours |
| is_fused_to(ring) | True if fused to another ring number |
| largest_number_of_bonds_shared_with_another_ring() | for flat ring systems, this will be 1 |
| strongly_fused_ring_neighbours() | Rings sharing more than 1 bond |
| contains_bond(a1, a2) | True if Ring contains these adjacent atoms |
| is_aromatic() | True if ring is aromatic |
| \__contains__ | Is atom included |
| \__len__ | Size |

To count the number of isolated (not fused) 5 membered aromatic rings
```
  result = 0
  m.compute_aromaticity_if_needed()

  for ring in m.rings():
    if len(ring) != 5:
      continue
    if ring.is_fused():
      continue;
    if ring.is_aromatic():
      result += 1
```

Counting the number of pyrrole type nitrogens, in isolated rings.
```
  result = 0
  for ring in m.rings():
    if len(ring) != 5:
      continue
    if ring.is_fused():
      continue;
    if not ring.is_aromatic():
      continue
    for atom in ring:
      if m.atomic_number(atom) == 7 && m.hcount(atom) == 1:
        result += 1
```
However it is unclear whether this would be more/or less efficient than the 
equivalent.

```
  result = 0
  for i, atom in enumerate(m):
    if atom.atomic_number() != 7:
      continue
    if not m.is_aromatic(i):
      continue
    if atom.hcount() == 0:
      continue
    if not m.in_ring_of_given_size(i, 5):
      continue
    if m.fused_system_size(i) == 1:
      result += 1
```

## Chemical Standardisation
Any work with molecules should ensure that molecules are represented in a consistent
manner. For example, are all the acids in charged or neutral forms? How are the nitro
groups represented? Etc...

Trying to formulate substruecture queries that can accommodate these variations is
challenging, and inefficient. LillyMol has a module that enforces consistent
molecular representations.

```
from lillymol_std import *
standardise = Chemical_Standardisation()
standardise.activate_all()
for mol in reader:
  standardise.process(m)
  # m is now in a consistent form for subsequent processing.
```

You may, or may not like how the molecules are changed, but they are all changed to
be consistent.

Within the C++ there are options for turning on just specific transformations, and
for transforming certain forms from LillyMol standard forms back to other forms;
transforming `N(=O)=O` to `[N+](=O)-[O-]` for example.

## Substructure Searching
LillyMol supports a rich set of substructure query capabilities. All invove a
`Substructure_Query` object that can be instantiate from

* smarts
* Molecule
* textproto file
* MDL query file

The current python implementation focuses on smarts and textproto forms.

Enable substructure searching via
```
from lillymol_query import *
```
Instantiate a new query for a para substituted methoxy group via the
one-line helper
```
query = QueryFromSmarts('[CH3]-[OD2]-c:c:c:[cD3]')
if query is None:
  logging.error('Invalid smarts')
```

This will generally be the most convenient way to build a query from SMARTS.
It returns a `SubstructureQuery` object on success and `None` if the SMARTS
cannot be parsed.

The equivalent two-step form is also available, which can be useful if you want
to allocate the object first and configure it explicitly
```
query = SubstructureQuery()
if not query.build_from_smarts('[CH3]-[OD2]-c:c:c:[cD3]'):
  logging.error('Invalid smarts')
```

To read a query from a textproto query specification
```
query = SubstructureQuery()
if not query.read_proto('/path/to/file.textproto'):
  logging.error('Cannot read query file %s...')
```

To perform a substructure search, not recording anything about the
atoms matched
```
m = MolFromSmiles('C(=N)(C1=C(O)C(=C(O)C=C1O)OC)CC1=CC=C(O)C=C1 CHEMBL503634')
query.substructure_search(m)
```
The number of matches will be returned. In this case it will frequently be 2 since
the query will match two times in a benzene like ring.

For a simple boolean test, use Python's `in` operator with the query on the
left and the molecule on the right
```
if query in m:
  print('matched')
```
This is equivalent to asking whether `query.substructure_search(m)` finds at
least one match, but it does not return the match count or matched atoms.

To get the matched atoms directly as a list of `Set_of_Atoms` objects, use:
```python
matches = query.substructure_search_matches(m)
```

For example if you wanted to place an isotope on each set of matched atoms that might look like
```
for match in query.substructure_search_matches(m):
  m.remove_isotopes()
  m.set_isotopes(match, 1)
  print(m)
```
Omit the `remove_isotopes` call to overwrite the new isotopes onto whatever isotopes might have already been there.

On the other hand if you need to know which matched atom is which, that might look like
```
for match in query.substructure_search_matches(m):
  m.remove_isotopes()
  for ndx, atom in enumerate(match):
    m.set_isotope(atom, ndx + 1)
  print(m)
```
Note that the isotope placed is incremented by 1 since isotope 0 does not mean anything.

The Substructure_Query class has a wide variety of options that control the matching. Those
are described in the `trxn` usage document. Here they are just listed

* set_only_keep_matches_in_largest_fragment
* set_embeddings_do_not_overlap
* set_find_one_embedding_per_atom
* set_find_unique_embeddings_only
* set_max_matches_to_find
* set_perceive_symmetry_equivalent_matches
* set_min_atoms_to_match
* set_max_atoms_to_match
* max_query_atoms_matched_in_search

### Examining atoms around a substructure match

A substructure search can populate an explicit `SubstructureResults` object. It
contains one or more embeddings, each represented by a `Set_of_Atoms` containing
the molecule atom numbers matched by the query. These embeddings can be passed
directly to `atoms_by_radius`:

```python
m = MolFromSmiles("CCOC")
query = QueryFromSmarts("[OD2]")
sresults = SubstructureResults()

if query.substructure_search(m, sresults):
  for matched_atoms in sresults:
    shells = m.atoms_by_radius(matched_atoms, 2)
    print("matched", list(shells[0]))
    print("one bond away", list(shells[1]))
    print("two bonds away", list(shells[2]))
```

`atoms_by_radius` performs one graph traversal and returns a list indexed by exact
minimum bond distance. Element zero is the starting set, element one contains
atoms one bond away, and so on. Each atom occurs in only one shell, even when it
is reachable from several matched atoms. Empty outer shells are retained. This
operation can therefore be embedded in a larger loop over molecules and query
embeddings without recomputing the inner shells for each radius.

## GFP fingerprint files and similarity search

`lillymol_tools.GFPList` provides read/search access to LillyMol GFP
fingerprints. It can read precomputed `.gfp`/TDT fingerprint files, and it can
also generate the standard LillyMol GFP fingerprint directly from Python
`Molecule` objects.

### Reading an existing GFP file

```python
from lillymol_tools import GFPList

gfp = GFPList.from_file("collection.gfp")

print(len(gfp))
print(gfp.tags())
print(gfp.smiles(0), gfp.id(0))

print(gfp.distance(0, 1))
```

The fingerprint schema is discovered from the tags in the first TDT in the file.
All subsequent fingerprints are interpreted with that same schema. Tag order in
the TDT is not significant, but all fingerprints in a `GFPList` must contain the
same components. Common tags include `MPR<` for molecular properties, `FP...<`
for fixed binary fingerprints, `NC...<` for sparse non-colliding fingerprints,
and `FC...<` for fixed counted fingerprints.

### Generating standard GFP fingerprints

The standard LillyMol GFP fingerprint can be generated directly from molecules.
This currently creates the long-standing standard GFP components:

- `FPIW<` path fingerprint
- `FPMK<` MACCS keys
- `FPMK2<` level-2 MACCS keys
- `MPR<` molecular properties

Use a `GFPContext` when you want standalone query fingerprints:

```python
from lillymol import Molecule
from lillymol_tools import GFPContext

ctx = GFPContext.standard()

mol = Molecule()
mol.build_from_smiles("CCO ethanol")
fp = ctx.fingerprint(mol)

print(ctx.distance(fp, fp))
```

The same standard context can also be built explicitly from immutable generator
specifications. This is the extension point for adding more fingerprint
generators while keeping the context schema explicit.

```python
from lillymol_tools import GFP, GFPContext

ctx = GFPContext.from_specs([
    GFP.iw(),
    GFP.maccs(level2=True),
    GFP.mpr(),
])
```

The component tags are generated by C++ from the full specification. Duplicate
tags are rejected, because a `GFPContext` cannot contain two indistinguishable
fingerprint components. The order of the specs is not significant; components
are canonicalized before the context hash is computed.

`GFP.maccs(level2=False)` generates only `FPMK<`; the default generates both
`FPMK<` and `FPMK2<`.

`GFP.alogp(replicates=9)` adds a sparse counted alogp component. The tag
contains the replicate count, for example `NCALOGP9<`, so contexts built with
different replicate counts are not compatible. The replicated bits all carry the
same bucketized alogp value as their count, matching the command-line alogp
fingerprint convention.

`GFP.xlogp(replicates=9)` adds the analogous sparse counted xlogp component
with tags such as `NCXLOGP9<`. It uses the xlogp command-line fingerprint
bucketization convention.

`GFP.tpsa(replicates=9)` adds a sparse counted topological polar surface area
component with tags such as `NCTPSA9<`. It preserves the command-line `psafp`
convention: TPSA is divided by 10, rounded to an integer bucket, and stored as
the count on each replicated bit, with a minimum count of 1.

`GFP.formula()` adds a fixed-size counted molecular formula fingerprint with
tag `FCFML<`. It uses the formula fingerprint log scaling defaults
recommended by the command-line tool, preserving the nonlinear count
scaling that emphasizes low element-count differences and compresses
differences among more abundant elements. The change going from one Oxygen
atom to two Oxygen atoms is much more important than going from 15 to 16
Carbon atoms.

`GFP.cats(max_path_length=10, include_hydrophobic_pairs=True)` adds the sparse
counted CATS pharmacophore-pair fingerprint. The tag is `NCCATS10<` when
hydrophobic-hydrophobic pairs are included and `NCCATSP10<` when they are
suppressed. The `P` suffix preserves the long-standing command-line convention
where the `jwcats -p` option suppresses hydrophobic-hydrophobic pairs. CATS
initialization requires `LILLYMOL_HOME` so the default charge and
donor/acceptor assigners can be loaded.

`GFP.ec(radius=3, atom_type="UST:Z")` adds a sparse counted extended
connectivity fingerprint. The default atom type is atomic number (`UST:Z`). The
radius and atom-type specification are encoded in the tag, for example
`NCEC3USTZ<` or `NCEC3USTAY<`. Colons in atom-type specifications are omitted
from tags; other non-alphanumeric characters are rejected. Multiple EC
fingerprints with different atom types can be combined in one context.

`GFP.ring_substitution()` adds the ring-substitution fingerprint with tag
`NCRS<`. It uses the defaults used by `gfp_make.pl`: full substituent atom
typing and single-feature bits.

Use `GFPList.standard()` when you want to build a searchable collection from
Python molecules:

```python
from lillymol import Molecule
from lillymol_tools import GFPContext, GFPList

gfp = GFPList.standard()

for smiles in ["CC ethane", "CCC propane", "CCCC butane"]:
  mol = Molecule()
  mol.build_from_smiles(smiles)
  gfp.add(mol)

query = Molecule()
query.build_from_smiles("CCC query")
query_fp = GFPContext.standard().fingerprint(query)

hits = gfp.nearest_neighbours(query_fp, 3)
for hit in hits:
  print(hit.index, hit.distance, gfp.id(hit.index))
```

If the molecules are already available as a Python list, build the list in one
call. By default this does not store smiles/id metadata because the caller
already owns the molecules and can keep whatever metadata structure is most
convenient. Set `store_metadata=True` if the `GFPList` should retain smiles and
ids for later `smiles(index)` and `id(index)` calls.

```python
molecules = []
for smiles in ["CC ethane", "CCC propane", "CCCC butane"]:
  mol = Molecule()
  mol.build_from_smiles(smiles)
  molecules.append(mol)

gfp = GFPList.standard_from_molecules(molecules)

# Later additions can also be batched. This also defaults to not storing
# smiles/id metadata.
gfp.add_molecules(more_molecules)
```

A `GFPList` either stores metadata for every entry or for none. Mixing the two
modes is rejected. Lists read from `.gfp` files and lists built with
`GFPList.add(mol)` store metadata. Batch construction defaults to no metadata.

Fingerprint generation follows the standard GFP preprocessing used by the
existing C++ GFP server path. The input Python molecule is copied before that
preprocessing, so calling `fingerprint()`, `GFPList.add()`, or
`GFPList.add_molecules()` does not alter the Python `Molecule` object. If
molecules have already been standardised and preprocessed, this step can be
skipped:

```python
ctx = GFPContext.standard(preprocess=False)
gfp = GFPList.standard(preprocess=False)
```

### Context compatibility

A `GFPFingerprint` is only meaningful with the `GFPContext` that generated it,
or with another context that has the same fingerprint schema. The Python binding
checks this compatibility before list/query distance calculations. Combining a
fingerprint generated by one context with a list using a different schema raises
`ValueError`.

This means a query fingerprint used against a `GFPList.standard()` collection
should be generated with `GFPContext.standard()` using the same settings. A
fingerprint generated with the standard context is not compatible with a list
read from a file that contains a different set of fingerprint tags.

### Nearest neighbours

Nearest-neighbour searches return `GFPNearestNeighbour` objects. The `index` is
the row number in the `GFPList`; use `smiles(index)` and `id(index)` to retrieve
the original molecule metadata.

```python
hits = gfp.nearest_neighbours(0, 10)
for hit in hits:
  print(hit.index, hit.distance, gfp.id(hit.index))

close = gfp.nearest_neighbours_within_distance(0, 0.25)
for hit in close:
  print(hit.index, hit.distance, gfp.smiles(hit.index), gfp.id(hit.index))
```

`nearest_neighbours(query, k)` returns up to `k` neighbours sorted by increasing
distance. `nearest_neighbours_within_distance(query, max_distance)` returns all
neighbours within the threshold, also sorted by increasing distance. When
`query` is a row number, the query fingerprint itself is not included in either
result. When `query` is a standalone `GFPFingerprint`, all list entries are
eligible neighbours.

### Component selection and weights

Fingerprint components can be selected and weighted by tag name.

```python
gfp.use_only(["FPIW<", "FPMK<"])
gfp.set_weight("FPMK<", 0.25)
gfp.use_all()
```

Weights affect subsequent distance and nearest-neighbour calculations. Changing
weights does not change the context compatibility hash, since the underlying
fingerprint schema is unchanged.

### Error handling

The Python API is designed to fail loudly in notebooks:

- unreadable GFP files raise `RuntimeError`
- trying to generate a fingerprint from a context that cannot generate
  fingerprints raises `RuntimeError`
- invalid row indices raise `IndexError`
- incompatible fingerprint contexts raise `ValueError`
- negative distance thresholds raise `ValueError`
- unknown component tags in `use_only` or `set_weight` raise `RuntimeError`

Repeated calls to `distance(i, j)` perform repeated distance computations. For
large all-against-all or repeated-nearest-neighbour workflows, prefer the
nearest-neighbour methods. A cached nearest-neighbour interface is planned for a
future version.

## Reactions
Enable reactions via
```
from lillymol_query import *
from lillymol_reaction import *
```
The `query` module must be imported first.

Reactions can be build from

* textproto reaction file
* smirks
* MDL reaction file.

```
  rxn = Reaction()
  if not rxn.read('/path/to/rxn.textproto'):
    logging.error('Cannot read reaction %s...
```
or
```
  rxn = Reaction()
  if not rxn.construct_from_smirks(smirks):
    logging.error('Cannot interpret smirks %s...
```
If the reaction is a simple form that has either no sidechains,
or all sidechains have a single, already specified, reagent, then
the `in_place_transformations` method can be used.
```
  rxn.in_place_transformations(m)
```
will perform a reaction to however many substructure matches
there are in `m`. This may, or may not be what you want. Since
a reaction inherits from a Substructure_Query object, there are
methods available for modifying the matching.

For more control over multiple matches, something like this may help
```
  rxn = Reaction()
  rxn.read('/path/to/reaction.textproto')
  matches = rxn.substructure_search_matches(m)
  # stop here if zero matches...
  for match in matches:
    product = rxn.perform_reaction(m, match)
    product.set_name(m.name() + ' ' + rxn.name())
    print(product)
```
or
```
  [product = rxn.perform_reaction(m, match) for match in matches]
```

### Multiple Reagents
The reaction object was designed to be able to rapidly enumerate large combinatorial
libraries. For this reason, it stores precomputed sets of sidechain reagents, which can
be rapidly added to a new scaffold. We cycle through these sidechains via an
iterator class. This workflow looks like 

1. Instantiate Reaction
2. Add reagents to the reaction
3. Process scaffols, generating multiple products for each scaffold.

In python, for a reaction with a single sidechain, processing a set of
molecules might look like
```
  rxn = Reaction()
  rxn.read('/path/to/reaction.textproto')
  rxn.add_sidechain_reagents(0, '/path/to/r2.smi', FileType.SMI)

  matches = rxn.substructure_search_matches(mol)
  if not matches:
    logging.error("No matches to %s", mol.name())

  iter = ReactionIterator(rxn)

  for scaffold in reader:
    iter.reset()
    while iter.active():
      for match in matches:
        product = rxn.perform_reaction(mol, match, iter)
        # do something with product

      iter++
```
The loop involving 'iter.active()' will loop over the reagents stored
in the reaction. For each such reagent, a product will be formed for
each embedding of the query in 'm', generating `number_reagents * number_matches`
products. Some molecules may have differing numbers of matches...

## Formal Charge Assignment
LillyMol contains a set of formal charge rules primarily developed by Robert F Bruns
at Lilly during the 1990's. These rules have proven robust over the years, with strong
concordance with other rule sets, with seemingly better results in many cases. From
the command line the script assign_formal_charges.sh assigns formal charges to a file
of molecules.

In order to access this functionality from python, a ChargeAssigner object must
be instantiated.
```
chg = ChargeAssigner()
for mol in mols:
  chg.process(mol)
```
If needed, the return code from `chg.process` is the number of formal charges assigned
to the molecule.

## Donor Acceptor Assigment.
LillyMol contains a set of donor acceptor rules primarily developed by Robert F Bruns
at Lilly during the 1990's. Acknowledge that definitions of donors and acceptors is
very complex with instances of various weak forms being found in various circumstances.
The rules apply mostly to strong donors and acceptors. One possibly novel concept is
the idea of a dual mode atom - something that can function either as a donor or an
acceptor. An OH atom is the most common example.

Results are returned by placing isotopic labels on the atoms. By convention isotope
1 is applied for acceptors, isotope 3 for donors, and isotope 2 to atoms which
can function as either.

A typical use might be
```
donor_acceptor = DonorAcceptor()
for mol in mols:
  donor_acceptor.process(mol)
```
The `process` method returns the number of assignments made. Note that if an atom
is both acceptor and donor, it will count twice in the return code value.

## Speed Comparisons
One speed comparison is described in [tsubstructure](/docs/Workflows/substructure_comparison.md).
That shows excellent speed from LillyMol python, but in that case there was not much
actually being done in python.

A more meaningful test was something to detect 4-pyridol groups and transform them to
4-pyridone types. Running across all of Chembl, this took about 6.5 minutes. Having
found an algorithm that worked, that was translated to C++ and runs in 1 minute 15
seconds. 

So, in cases where most of the calculation is being done with python, speed
may be significantly diminished.

## Miscellaneous functions

### count_atoms_in_smiles
Many times the only need to instantiate a molecule is to get an atom count.
This function is text only, and makes, what is usually a very good, count of the
number of heavy atoms in a smiles.

### set_auto_create_new_elements
Allow creation of arbitrary elements
```
[Th][Eq][U]IC[K][Br]O[W]NFO[Xj][U][MP]SO[Ve][R][Th][E][La][Zy][D]O[G]
```
to be a valid smiles. Substructure search this with
```
[#{Th}][#{Eq}]...
``` 

### set_atomic_symbols_can_have_arbitrary_length
Allow molecules such as
```
[Ala][Gly][Arg][Ser]
```
again substructure searching is with `[#{Ala}D1]~[#{Gly}D2]...`.

### interpret_D_as_deuterium, interpret_T_as_deuterium
Elements `D` and `T` are interpreted as Hydrogen isotopes.

