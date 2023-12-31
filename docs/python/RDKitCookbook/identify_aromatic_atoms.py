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
