"""Time the formula fingerprint as a filter ahead of a more precise, slower, second stage.

The formula fingerprint is not selective, so many molecules that do not contain a query
pass it. Its advantage is speed: it is much faster to compare than the second stage, so
if it removes enough molecules it saves time, and the second stage need only be applied
to the molecules that pass.

The second stage is, by default, a substructure search. With --second-stage it can instead
be a richer fingerprint, made from the smiles for each molecule, and compared with the same
Tversky similarity: the path fingerprint (GFP.iw()) or the extended connectivity
fingerprint (GFP.ec()).

This script measures that on a collection of molecules. It
  1. reads the formula fingerprints of the collection and times the Tversky search,
  2. times the second stage on a random sample of the molecules,
  3. checks on that sample that the filter did not reject any molecule that the second
     stage finds, and
  4. estimates the time to search the whole collection with and without the filter.

With numpy the filter uses gfp.indices_within_distance, which returns the row numbers as an
array. Without it the filter uses gfp.nearest_neighbours_within_distance, which makes an object
for each molecule that passes, and takes much longer when many pass.

Usage: formula_prefilter_timing.py [FILE.gfp] [--query SMILES] [--sample N]
                                   [--second-stage substructure|path|ec]
                                   [--second-stage-distance D]

FILE.gfp contains the fingerprints and smiles of the collection, without the counts of
atoms bearing hydrogens, which would prevent the filter finding molecules that contain
the query. Make it with
  formula_fingerprint -H collection.smi > collection.formula_nh.gfp
or, if the molecules are not in a smiles file,
  gfp_make.sh -FORMULA,nohydrogen collection.smi > collection.formula_nh.gfp
With no file a small bundled collection is used, which is far too small for the times to
mean anything.

The query is used as a SMILES to make the fingerprints and, for the substructure second
stage, as a SMARTS, so it must be valid as both.

With a fingerprint as the second stage, a molecule is found if Tversky(query=2, candidate=0)
gives a distance of at most --second-stage-distance (default 0), which means the molecule
has at least the features of the query. If the filter rejects some of those molecules that is
reported but is not an error, since two fingerprints can disagree. With a substructure search
a molecule that contains the query must never be rejected, so that is an error.
"""
import argparse
import random
import time
from pathlib import Path

try:
  import numpy as np
except ImportError:  # The filter then makes an object for each result, which is slower.
  np = None

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


def make_second_stage(kind, query, distance):
  """Returns a description, a word for what it finds, and a function of a molecule that
  returns true if the second stage finds it."""
  if kind == "substructure":
    substructure = SubstructureQuery()
    if not substructure.build_from_smarts(query):
      raise SystemExit(f"Invalid query '{query}'")
    return ("substructure search", "contain the query",
            lambda molecule: substructure.substructure_search(molecule) > 0)

  spec, name = (GFP.iw(), "path fingerprint") if kind == "path" else (GFP.ec(), "EC fingerprint")
  context = GFPContext.from_specs([spec])
  query_fp = context.fingerprint(MolFromSmiles(query + " query"))
  tversky = Tversky(query=2.0, candidate=0.0)
  limit = distance + 1.0e-6
  # The fingerprint of every molecule has to be made, which is most of the cost.
  return (f"{name} made, then Tversky(2, 0)", f"have a distance of at most {distance:g}",
          lambda molecule: context.distance(query_fp, context.fingerprint(molecule), tversky) <= limit)


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
                      help="molecules used to time the second stage (default 20000)")
  parser.add_argument("--second-stage", choices=["substructure", "path", "ec"],
                      default="substructure",
                      help="what is done to the molecules that pass the filter (default "
                           "substructure): a substructure search, or make a path or EC "
                           "fingerprint and compare it with Tversky(2, 0)")
  parser.add_argument("--second-stage-distance", type=float, default=0.0,
                      help="with a fingerprint second stage, the largest distance that is "
                           "found (default 0)")
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

  # The scan itself. Asking for the one nearest neighbour compares the query with every
  # molecule, but builds almost no list of results.
  _, scan_time = timed(lambda: gfp.nearest_neighbours(query_fp, 1, tversky=contains_query),
                       repeats=3)

  # The filter: everything with at least the counts of the query, as an array of the row
  # numbers. This makes no object for each molecule that passes, which would be most of the
  # time when many pass (nearest_neighbours_within_distance does that).
  if np is not None:
    rows, filter_time = timed(
        lambda: gfp.indices_within_distance(query_fp, 1.0e-6, tversky=contains_query),
        repeats=3)
    passed = np.zeros(n, dtype=bool)
    passed[rows] = True
  else:
    hits, filter_time = timed(
        lambda: gfp.nearest_neighbours_within_distance(query_fp, 1.0e-6, tversky=contains_query),
        repeats=3)
    rows = [hit.index for hit in hits]
    passed = [False] * n
    for row in rows:
      passed[row] = True
  print(f"\nFilter, Tversky(query=2, candidate=0) at distance 0, query {args.query}")
  print(f"  scan only        {format_time(scan_time)}  ({scan_time / n * 1.0e9:.0f} ns per molecule)")
  print(f"  scan and rows    {format_time(filter_time)}  ({filter_time / n * 1.0e9:.0f} ns per molecule)")
  print(f"  passing          {len(rows)} of {n} molecules ({100.0 * len(rows) / n:.1f}%)")

  # The second stage, timed on a sample, since it is much slower per molecule.
  description, finds, found_by_second_stage = make_second_stage(
      args.second_stage, args.query, args.second_stage_distance)
  sample = random.Random(1).sample(range(n), min(args.sample, n))

  # The smiles has to be parsed, as it would be in any real use. If the molecules are
  # already in memory only the time for the second stage applies.
  molecules, parse_time = timed(lambda: [MolFromSmiles(gfp.smiles(i)) for i in sample])
  found, stage_time = timed(
      lambda: [i for i, m in zip(sample, molecules) if found_by_second_stage(m)])
  per_molecule = (parse_time + stage_time) / len(sample)
  print(f"\nSecond stage, {description}, of {len(sample)} molecules")
  print(f"  parse the smiles   {format_time(parse_time)}  "
        f"({parse_time / len(sample) * 1.0e6:.1f} us per molecule)")
  print(f"  {args.second_stage:18s} {format_time(stage_time)}  "
        f"({stage_time / len(sample) * 1.0e6:.1f} us per molecule)")
  print(f"  {finds}  {len(found)} of {len(sample)}")

  rejected = [i for i in found if not passed[i]]
  print(f"  rejected by the filter  {len(rejected)} of those found by the second stage")

  # Estimates for the whole collection. The filter time is the one that includes the
  # row numbers, which is what is needed to know which molecules to search.
  without_filter = per_molecule * n
  with_filter = filter_time + per_molecule * len(rows)
  print(f"\nEstimated time to search all {n} molecules")
  print(f"  second stage only            {format_time(without_filter)}")
  print(f"  filter, then second stage    {format_time(with_filter)}  "
        f"({without_filter / with_filter:.1f} times faster)")
  print(f"The filter cannot make it more than {n / max(len(rows), 1):.1f} times faster, "
        f"however slow the second stage is, since the molecules that pass still have to be "
        f"processed.")
  print("These are estimates from the sample. The filter is not selective, so how much "
        "it saves depends on the query.")

  if bundled:
    assert len(rows) == 21  # The same as formula_tversky.py for a phenol query.
    if args.second_stage == "substructure":
      assert not rejected
  if rejected and args.second_stage == "substructure":
    raise SystemExit(f"{len(rejected)} molecules that contain the query were rejected, "
                     f"was the file made without the hydrogen counts?")

if __name__ == "__main__":
  main()
