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


# A reaction involving multiple sidechains as substituents.

rxn_text = '''
scaffold {
  id: 0
  smarts: "[#0]c"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: "[CD<2]"
  match_conditions {
    ignore_symmetry_related_matches: true
  }
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
'''
rxn = Reaction()
rxn.construct_from_textproto(rxn_text)

scaffold = MolFromSmiles("*c1ccccc1")
sidechains = [MolFromSmiles(s) for s in ["C", "CC", "CCC", "CC1CC1", "CCCC"]]

smc = SidechainMatchConditions()
for sidechain in sidechains:
  rxn.add_sidechain_reagent(0, sidechain, smc)

iter = ReactionIterator(rxn)
sresults = SubstructureResults();
rxn.substructure_search(scaffold, sresults)
for embedding in sresults:
  while iter.active():
    product = rxn.perform_reaction(scaffold, embedding, iter)
    print(product.aromatic_smiles())
    iter.increment()
