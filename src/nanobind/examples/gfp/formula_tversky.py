"""Tversky searching with the molecular formula fingerprint.

The formula fingerprint is a fixed size counted fingerprint. Its slots count
carbon, nitrogen, oxygen and sulphur (aromatic and aliphatic separately), the
other common elements, the number of atoms that carry a hydrogen (on
heteroatoms, aromatic carbons and other carbons), and ring atoms. The counts are
scaled nonlinearly so differences between small counts matter most.

Tversky similarity is asymmetric, so it can ask "which molecules have at least
as much of everything as my query?" or "which molecules have nothing that my
query lacks?". Those are the two Tversky searches below.

Usage: formula_tversky.py [haystack.smi]

With no argument the bundled molecules.smi is searched and the results are
checked. Any other file is searched and the results are just printed.
"""
import sys
from pathlib import Path

from lillymol import GFP, GFPContext, GFPList, MolFromSmiles, Tversky, slurp


def search(gfp, query_fp, tversky, k):
  """Returns [(id, distance)] for the k nearest neighbours of query_fp."""
  hits = gfp.nearest_neighbours(query_fp, k, tversky=tversky)
  return [(gfp.id(hit.index), hit.distance) for hit in hits]


def main(argv):
  default_file = Path(__file__).with_name("molecules.smi")
  smiles_file = argv[1] if len(argv) > 1 else str(default_file)

  molecules = slurp(smiles_file)

  # A context that generates only the formula fingerprint (tag FCFML<).
  context = GFPContext.from_specs([GFP.formula()])

  # Fingerprint everything. store_metadata keeps the smiles and ids.
  gfp = GFPList(context)
  gfp.add_molecules(molecules, store_metadata=True)
  print(f"{len(gfp)} molecules fingerprinted from {smiles_file}")

  # A made up query, here paracetamol, C8H9NO2.
  query_fp = context.fingerprint(MolFromSmiles("CC(=O)Nc1ccc(O)cc1 query"))

  k = 5

  # Tanimoto: no Tversky object.
  tanimoto = search(gfp, query_fp, None, k)

  # Tversky(query=a, candidate=b). `query` weights what is only in the query,
  # `candidate` weights what is only in the molecule being compared. Giving all
  # the weight to the query means anything the candidate has in addition is free,
  # so molecules with at least as much of everything as the query are at distance 0.
  contains_query = search(gfp, query_fp, Tversky(query=2.0, candidate=0.0), k)

  # The other way round: molecules with nothing that the query lacks.
  contained_in_query = search(gfp, query_fp, Tversky(query=0.0, candidate=2.0), k)

  for title, hits in (("Tanimoto", tanimoto),
                      ("Tversky query=2 candidate=0, candidates with at least what the query has",
                       contains_query),
                      ("Tversky query=0 candidate=2, candidates with nothing the query lacks",
                       contained_in_query)):
    print(f"\n{title}")
    for identifier, distance in hits:
      print(f"  {distance:6.3f}  {identifier}")

  if len(argv) == 1:
    # Check the bundled results.
    zero = lambda hits: {name for name, distance in hits if distance < 1.0e-6}
    assert zero(tanimoto) == {"paracetamol"}
    assert zero(contains_query) == {"paracetamol", "propionyl_paracetamol",
                                    "chloroacetyl_paracetamol"}
    assert zero(contained_in_query) == {"paracetamol", "aminophenol", "acetamide"}


if __name__ == "__main__":
  main(sys.argv)
