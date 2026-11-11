"""Time the formula fingerprint as a filter ahead of substructure searching.

The formula fingerprint is not selective, so many molecules that do not contain a query
pass it. Its advantage is speed: it is much faster to compare than to do a substructure
search, so if it removes enough molecules it saves time, and the substructure search
need only be applied to the molecules that pass.

This script measures that on a collection of molecules. It
  1. reads the formula fingerprints of the collection and times the Tversky search,
  2. times the precise (substructure) search on a random sample of the molecules,
  3. checks on that sample that the filter did not reject any molecule that contains the
     query, and
  4. estimates the time to search the whole collection with and without the filter.

Usage: formula_prefilter_timing.py [FILE.gfp] [--query SMILES] [--sample N]

FILE.gfp contains the fingerprints and smiles of the collection, without the counts of
atoms bearing hydrogens, which would prevent the filter finding molecules that contain
the query. Make it with
  formula_fingerprint -H collection.smi > collection.formula_nh.gfp
or, if the molecules are not in a smiles file,
  gfp_make.sh -FORMULA,nohydrogen collection.smi > collection.formula_nh.gfp
With no file a small bundled collection is used, which is far too small for the times to
mean anything.

The query is used as a SMILES to make the fingerprint and as a SMARTS for the precise
search, so it must be valid as both.
"""
import argparse
import random
import time
from pathlib import Path

from lillymol import (GFP, GFPContext, GFPList, MolFromSmiles, SubstructureQuery,
                      Tversky)

TAG = "FCFMLNH<"


def timed(function, repeats=1):
  """Returns the result of `function` and the shortest time, in seconds, of `repeats` calls."""
  best = None
  for _ in range(repeats):
    start = time.perf_counter()
    result = function()
    elapsed = time.perf_counter() - start
    best = elapsed if best is None else min(best, elapsed)
  return result, best


def format_time(seconds):
  if seconds < 1.0e-3:
    return f"{seconds * 1.0e6:.0f} us"
  if seconds < 1.0:
    return f"{seconds * 1.0e3:.1f} ms"
  if seconds < 120.0:
    return f"{seconds:.2f} s"
  return f"{seconds / 60.0:.1f} min"


def main():
  parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
  parser.add_argument("file", nargs="?", help="formula fingerprints, see above")
  parser.add_argument("--query", default="c1ccccc1O", help="query SMILES (default: phenol)")
  parser.add_argument("--sample", type=int, default=20000,
                      help="molecules used to time the precise search (default 20000)")
  args = parser.parse_args()

  bundled = args.file is None
  fname = args.file or str(Path(__file__).with_name("molecules.formula_nh.gfp"))

  gfp, read_time = timed(lambda: GFPList.from_file(fname))
  n = len(gfp)
  print(f"{n} molecules read from {fname} in {format_time(read_time)}")
  if gfp.tags() != [TAG]:
    raise SystemExit(f"The file has fingerprints {gfp.tags()} but {TAG} is needed, make "
                     f"it with formula_fingerprint -H. Without the hydrogen counts a "
                     f"molecule that contains the query has at least its counts.")
  if not gfp.metadata_stored():
    raise SystemExit("The file must contain the smiles of the molecules")

  # The query. The fingerprint must be made with the same settings as the file.
  context = GFPContext.from_specs([GFP.formula(hydrogens=False)])
  query_fp = context.fingerprint(MolFromSmiles(args.query + " query"))
  contains_query = Tversky(query=2.0, candidate=0.0)

  # The filter: everything with at least the counts of the query. The time includes
  # making the list of results.
  hits, filter_time = timed(
      lambda: gfp.nearest_neighbours_within_distance(query_fp, 1.0e-6, tversky=contains_query),
      repeats=3)
  passed = {hit.index for hit in hits}
  print(f"\nFilter, Tversky(query=2, candidate=0) at distance 0, query {args.query}")
  print(f"  time      {format_time(filter_time)}  ({filter_time / n * 1.0e9:.0f} ns per molecule)")
  print(f"  passing   {len(passed)} of {n} molecules ({100.0 * len(passed) / n:.1f}%)")

  # The precise search, timed on a sample, since it is much slower per molecule.
  substructure = SubstructureQuery()
  if not substructure.build_from_smarts(args.query):
    raise SystemExit(f"Invalid query '{args.query}'")
  sample = random.Random(1).sample(range(n), min(args.sample, n))

  # The smiles has to be parsed, as it would be in any real use. If the molecules are
  # already in memory only the time for the substructure search applies.
  molecules, parse_time = timed(lambda: [MolFromSmiles(gfp.smiles(i)) for i in sample])
  contain, search_time = timed(
      lambda: [i for i, m in zip(sample, molecules) if substructure.substructure_search(m) > 0])
  per_molecule = (parse_time + search_time) / len(sample)
  print(f"\nPrecise search of {len(sample)} molecules")
  print(f"  parse the smiles   {format_time(parse_time)}  "
        f"({parse_time / len(sample) * 1.0e6:.1f} us per molecule)")
  print(f"  substructure search {format_time(search_time)}  "
        f"({search_time / len(sample) * 1.0e6:.1f} us per molecule)")
  print(f"  contain the query  {len(contain)} of {len(sample)}")

  rejected = [i for i in contain if i not in passed]
  print(f"  rejected by the filter  {len(rejected)} of those that contain the query")

  # Estimates for the whole collection.
  without_filter = per_molecule * n
  with_filter = filter_time + per_molecule * len(passed)
  print(f"\nEstimated time to search all {n} molecules")
  print(f"  precise search only          {format_time(without_filter)}")
  print(f"  filter, then precise search  {format_time(with_filter)}  "
        f"({without_filter / with_filter:.1f} times faster)")
  print("These are estimates from the sample. The filter is not selective, so how much "
        "it saves depends on the query.")

  if bundled:
    assert not rejected
    assert len(passed) == 21  # The same as formula_tversky.py for a phenol query.
  if rejected:
    raise SystemExit(f"{len(rejected)} molecules that contain the query were rejected, "
                     f"was the file made without the hydrogen counts?")


if __name__ == "__main__":
  main()
