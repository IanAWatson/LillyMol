from lillymol import *


mol = MolFromSmiles("c1cccc2c1CCCC2")
mol.compute_aromaticity_if_needed()
for ring in mol.rings():
  print(ring)
  print(ring.is_aromatic())
