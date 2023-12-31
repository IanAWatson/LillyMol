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
The case used in the RDKit Cookbook is conceptually different from
how such things would typically be processed in LillyMol. This example
needs to be converted to atom numbers, because if we remove a bond number,
that may invalidate any exisitng bond numbers.

So this example is re-cast to use atom numbers. There is no library
function to perform this operation, although clearly it would be
straightforward.
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

When molecules
are broken, that is typically the result of a substructure search and
things can be done by atom numbers, not bond numbers. 
```
