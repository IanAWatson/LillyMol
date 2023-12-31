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

