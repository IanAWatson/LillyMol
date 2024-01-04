# RDKit Cookbook
The RDKit Cookbook contains useful and interesting recipes for performing common
tasks with the RDKit python interface. Many relate to structure depictions,
but some of those that are more generic are repeated here.

# Bonds and Bonding

## Hybridization Type and Count
[hybridization_type_and_count](hybridization_type_and_count.py)
This example extends the RDKit example to count the number of atoms with
each hybridization state. We also show how the result can be obtained
by interrogation either the molecule or the individual atoms.
```
from collections import defaultdict

from lillymol import *

def hybridization(unsaturation:int)->str:
  """Return a string for the Hybridziation of an atom with unsaturation
    Args:
      unsaturation: the difference btw nbonds and ncon for an atom
    Returns:
      String, SP, SP2, SP3
  """
  return ["SP3", "SP2", "SP"][unsaturation]

m = MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C#C")

# Accumulate the number of atoms with each hybridization type.
total = defaultdict(int)

for i in range(0, len(m)):
  total[hybridization(m.unsaturation(i))] += 1

for (hyb,count) in total.items():
  print(f"{hyb} {count}")

# The atom object also has an unsaturation() method so we can access
# values via that means

total = defaultdict(int)
for atom in m:
  total[hybridization(atom.unsaturation())] += 1

for (hyb,count) in total.items():
  print(f"{hyb} {count}")

```

## Count Ring Systems
[count_ring_systems](count_ring_systems.py)
This example shows how LillyMol provides some higher lavel abstractions of
ring systems, which can lead to considerable simplification and efficiencies.
```
from lillymol import *

def get_ring_systems(mol:Molecule, include_spiro=False):
  # Return a list of the atoms in each ring system.
  if include_spiro:
    sysid = m.label_atoms_by_ring_system_including_spiro_fused()
  else:
    sysid = m.label_atoms_by_ring_system()

  result = []
  nsys = max(sysid)

  # non ring atoms get ring system id 0 so start at 1.
  for sys_num in range(1, nsys + 1):
    result.append([i for i,s in enumerate(sysid) if s == sys_num])

  return result

m = MolFromSmiles("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3")
rsys = get_ring_systems(m)
for r in rsys:
  print(set(r))
```

## Identify Aromatic Rings
[identify_aromatic_rings](identify_aromatic_rings.py)
Aromaticity of rings is stored as a property of the ring objects, and so the
aromaticity of a ring is a property of that ring.

There is
no real concept of iterating through the bond numbers that define a ring. In fact
there are very few cases where iterating through a list of bond numbers is used.
It is however fairly common to iterate through a ring while keeping track of
the previous atom, which does define a bond since the atoms in a Ring are stored
in order, and adjacent atoms are bonded - and the first and last items are bonded.
```
from lillymol import *


mol = MolFromSmiles("c1cccc2c1CCCC2")
mol.compute_aromaticity_if_needed()
for ring in mol.rings():
  print(ring)
```

## Identify Aromatic Atoms
[identify_aromatic_atoms](identify_aromatic_atoms.py)
The G smarts extension is unsaturation, nbonds() - ncon(), so a fully 
saturated atom matches `G0`, an olefin or aromatic carbon will match `G1`
and a cyano will match `G2`.
```
from lillymol import *
from lillymol_query import *

mol = MolFromSmiles("c1ccccc1C=CCC")
aromatic_carbon = SubstructureQuery()
aromatic_carbon.build_from_smarts("c")
sresults = SubstructureResults()
aromatic_carbon.substructure_search(mol, sresults)
print("Aromatic carbon atoms")
for embedding in sresults:
  print(embedding)

# The G smarts extension means unsaturation.
olefinic_carbon = SubstructureQuery()
olefinic_carbon.build_from_smarts("[CG1]")
olefinic_carbon.substructure_search(mol, sresults)
print("olefinic carbon")
for embedding in sresults:
  print(embedding)
```

## Create Fragments
In Lillymol the way to do this is with command line tool [dicer](/docs/Molecule_Tools/dicer.md),
which offers a lot of flexibility and performance, and addresses the recursive
fragment creation problem. In any situation where molecular fragmentation
was needed, `dicer` would be the tool to use.

The case used in the RDKit Cookbook is conceptually different from
how such things would typically be processed in LillyMol. This example
needs to be converted to atom numbers, because if we remove a bond by number,
that may invalidate any exisitng bond numbers - inside a Molecule the
bonds are stored as a vector of bonds, so removing any one can
invalidate previous atom numbers.

This code seems to recreate what is done.
```
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
```

A simple recusrive version of fragment breaking below generates 74 thousand
fragments from the starting molecule below.
```

def create_fragments(mol:Molecule):
  mol.compute_aromaticity_if_needed()
  print(mol.smiles())
  breakable_bonds = []

  for bond in mol.bonds():
    if bond.nrings():
      continue
    if not bond.is_single_bond():
      continue
    breakable_bonds.append(bond)

  for bond in breakable_bonds:
    mol.remove_bond_between_atoms(bond.a1(), bond.a2())
    for component in mol.create_components():
      create_fragments(component);
    mol.add_bond(bond.a1(), bond.a2(), SINGLE_BOND)

mol = MolFromSmiles("CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C([N-]2)C=CC(=C3)OC")
create_fragments(mol)
```
takes 0.7 seconds. The reason for the large number of fragments is that the
molecule contains 9 breakable bonds.

## Largest Fragment
Just as with RDKit there are several ways of doing this.
```
from lillymol import *

mol = MolFromSmiles("CCOC(=O)C(C)(C)OC1=CC=C(C=C1)Cl.CO.C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N")
mol.reduce_to_largest_fragment_carefully()
print(mol.aromatic_smiles())

mol = MolFromSmiles("CCOC(=O)C(C)(C)OC1=CC=C(C=C1)Cl.CO.C1=CC(=CC=C1C(=O)N[C@@H](CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N")
components = mol.create_components()
largest = max(components, key=lambda m: m.natoms())
print(largest.aromatic_smiles())
```

## Sidechain Core Enumeration
Generally this can be done most conveniently at the command line with [trxn](/docs/Molecule_Tools/trxn.md).
Smirks support in LillyMol covers most simple cases, but there are limitations. In the case of this
reaction, in order to have the `*` atoms removed, they must be assigned atom labels.
```
from lillymol import *
from lillymol_query import *
from lillymol_reaction import *

core = MolFromSmiles("*c1c(C)cccc1O")
sidechain = MolFromSmiles("CN*")

rxn = Reaction()
rxn.construct_from_smirks("[c:1][#0:3].[#0:4][*:2]>>[*:1]-[*:2]")

products = rxn.perform_reaction(core, sidechain)
for product in products:
  print(product.unique_smiles())
```
Rather than smirks, the reaction can be constructed from a textproto reaction specification.
```
rxn_text = '''
scaffold {
  id: 0
  smarts: "[#0]c"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: "[#0]*"
  remove_atom: 0
  join {
    a1: 1
    a2: 1
    btype: SS_SINGLE_BOND
  }
}
'''
rxn = Reaction()
rxn.construct_from_textproto(rxn_text)

products = rxn.perform_reaction(core, sidechain)
for product in products:
  print(product.unique_smiles())
```
yields the same results. The textproto reaction specification is much more
verbose, but very precise and directive in what changes are performed. And
of course, query files can be re-used for reactions this way.

Reaction objects in LillyMol were originally designed to process combinatorial
libraries, so they naturally accommodate the idea of adding multiple
sidechains to a scaffold - library enumeration would typically involve
processing one scaffold at a time, adding all sidechains to it.
