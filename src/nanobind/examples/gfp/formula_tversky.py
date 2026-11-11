"""Approximate substructure searching with Tversky similarity and the formula fingerprint.

The formula fingerprint is a fixed size counted fingerprint: it counts atoms of each
kind (carbon, nitrogen, oxygen and sulphur, each aromatic or not, the other common
elements, and ring atoms). It is very fast to generate and to compare, but it knows
nothing about how atoms are connected, so it cannot say that a molecule contains a
query, only that it has at least as many atoms of each kind as the query.
That makes it a coarse first stage: apply it to a large collection to find the
molecules that might contain the query, then use substructure searching, or another
precise method, on just those.

Tversky similarity is asymmetric, so it can ask exactly that question:

  Tversky(query=2, candidate=0)  distance 0 if the candidate has at least every count
                                 of the query, however many more atoms it has.
  Tversky(query=0, candidate=2)  distance 0 if the candidate has nothing the query lacks.
  Tversky(query=1, candidate=1)  Tanimoto.

This script
  1. shows how the ranking of molecules changes as the Tversky weights change, and
  2. checks how well the fingerprint works as a filter, against real substructure
     searching, with and without the counts of atoms bearing hydrogens.

Usage: formula_tversky.py [--file molecules.smi] [--query SMILES]

The query is used as a SMILES to make the fingerprint and as a SMARTS for the precise
search, so it must be valid as both. With no arguments a bundled set of molecules and a
phenol query are used, and the results are checked.
"""
import argparse
from pathlib import Path

from lillymol import (GFP, GFPContext, GFPList, MolFromSmiles, SubstructureQuery,
                      Tversky, slurp)

# (query weight, candidate weight, label). Q is the query and C is the candidate.
SETTINGS = [
    (1.0, 1.0, "Tanimoto"),
    (1.5, 0.5, "mostly Q"),
    (1.9, 0.1, "nearly Q"),
    (2.0, 0.0, "contains Q"),
    (0.5, 1.5, "mostly C"),
    (0.0, 2.0, "inside Q"),
]

# Molecules to follow through the settings, if they are in the file.
PROBES = ["phenol", "o-cresol", "bisphenol_A", "benzyl_alcohol", "benzene", "cyclohexanol",
          "ethanol", "glycine"]


def build(smiles_file, hydrogens):
  """Returns the context, and a GFPList of the fingerprints of the molecules in the file."""
  context = GFPContext.from_specs([GFP.formula(hydrogens=hydrogens)])
  gfp = GFPList(context)
  gfp.add_molecules(slurp(smiles_file), store_metadata=True)
  return context, gfp


def distances(gfp, query_fp, query, candidate):
  """The distance from the query to every molecule, as {identifier: distance}."""
  tversky = Tversky(query=query, candidate=candidate)
  hits = gfp.nearest_neighbours(query_fp, len(gfp), tversky=tversky)
  return {gfp.id(hit.index): hit.distance for hit in hits}


def sweep(gfp, query_fp):
  """Show how the distances to some molecules change with the Tversky weights."""
  by_setting = [distances(gfp, query_fp, q, c) for q, c, _ in SETTINGS]

  names = [name for name in PROBES if name in by_setting[0]]
  if len(names) < len(PROBES):  # Not the bundled molecules, use the nearest by Tanimoto.
    names = sorted(by_setting[0], key=by_setting[0].get)[:len(PROBES)]

  print("Distance from the query, as the Tversky weights change")
  print("(query weight, candidate weight. A weight of 0 means that kind of difference is free.")
  print(" Q is the query and C is the candidate. 'contains Q' is the filter used below.)\n")
  print(f"{'':16s}" + "".join(f"{f'{q:g}, {c:g}':>12s}" for q, c, _ in SETTINGS))
  print(f"{'':16s}" + "".join(f"{label:>12s}" for _, _, label in SETTINGS))
  for name in names:
    print(f"{name:16s}" + "".join(f"{d[name]:12.3f}" for d in by_setting))
  print("\nThe fingerprint does not know how atoms are connected, so molecules with the same")
  print("atoms, such as o-cresol and benzyl_alcohol, cannot be told apart. Aromatic and")
  print("aliphatic atoms are counted separately, which is why benzene is 'inside' phenol")
  print("but ethanol, with its aliphatic carbons, is not.\n")
  return by_setting


def precise_matches(smiles_file, query):
  """The identifiers of the molecules that contain `query` as a substructure."""
  substructure = SubstructureQuery()
  if not substructure.build_from_smarts(query):
    raise SystemExit(f"Invalid query '{query}'")
  return {str(m.name()) for m in slurp(smiles_file) if substructure.substructure_search(m) > 0}


def filter_check(smiles_file, query, truth):
  """For each form of the fingerprint, what the filter finds, against `truth`."""
  print(f"Used as a filter, Tversky(query=2, candidate=0) at distance 0, "
        f"query {query}")
  print(f"{len(truth)} molecules contain the query\n")
  results = {}
  for hydrogens in (False, True):
    context, gfp = build(smiles_file, hydrogens)
    query_fp = context.fingerprint(MolFromSmiles(query + " query"))
    # Everything that has at least the counts of the query.
    hits = gfp.nearest_neighbours_within_distance(
        query_fp, 1.0e-6, tversky=Tversky(query=2.0, candidate=0.0))
    passed = {gfp.id(hit.index) for hit in hits}
    results[hydrogens] = passed

    print(f"GFP.formula(hydrogens={hydrogens}), tag {context.tags()[0]}")
    print(f"  passes the filter  {len(passed):4d} of {len(gfp)} molecules")
    print(f"  contain the query  {len(truth & passed):4d} of {len(truth)} "
          f"({100.0 * len(truth & passed) / max(len(truth), 1):.0f}% found)")
    if truth - passed:
      print(f"  missed             {', '.join(sorted(truth - passed))}")
    if passed - truth:
      print(f"  do not contain it  {', '.join(sorted(passed - truth))}")
    print()
  return results


def main():
  parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
  parser.add_argument("--file", help="SMILES file to search (default: the bundled molecules)")
  parser.add_argument("--query", default="c1ccccc1O", help="query SMILES (default: phenol)")
  args = parser.parse_args()

  bundled = args.file is None
  smiles_file = args.file or str(Path(__file__).with_name("molecules.smi"))

  # Without the hydrogen counts a molecule that contains the query has at least its counts.
  context, gfp = build(smiles_file, hydrogens=False)
  print(f"{len(gfp)} molecules from {smiles_file}, query {args.query}\n")
  query_fp = context.fingerprint(MolFromSmiles(args.query + " query"))

  by_setting = sweep(gfp, query_fp)

  truth = precise_matches(smiles_file, args.query)
  results = filter_check(smiles_file, args.query, truth)

  if bundled and args.query == "c1ccccc1O":
    zero = {name for name, d in by_setting[3].items() if d < 1.0e-6}
    # What the filter finds is the same as what Tversky(2, 0) gives at distance zero.
    assert zero == results[False]
    # Without the hydrogen counts nothing that contains the query is lost.
    assert truth <= results[False]
    # It is not selective, some molecules pass that do not contain the query.
    assert "benzyl_alcohol" in results[False] - truth
    assert "cyclohexanol" not in results[False]
    # With the hydrogen counts most of the molecules that contain the query are lost.
    assert len(truth & results[True]) < len(truth) / 2
    assert {"o-cresol", "paracetamol", "aspirin"} <= truth - results[True]
    # At weights 1, 1 it is Tanimoto: only the molecule itself has distance 0.
    assert {name for name, d in by_setting[0].items() if d < 1.0e-6} == {"phenol"}
    # At weights 0, 2 only molecules with nothing the query lacks have distance 0.
    # Aromatic and aliphatic carbon are counted separately, so ethanol is not inside.
    inside = {name for name, d in by_setting[5].items() if d < 1.0e-6}
    assert {"phenol", "benzene"} <= inside
    assert "o-cresol" not in inside and "ethanol" not in inside
    # Molecules with the same atoms cannot be told apart, whatever the weights.
    for d in by_setting:
      assert abs(d["o-cresol"] - d["benzyl_alcohol"]) < 1.0e-6


if __name__ == "__main__":
  main()
