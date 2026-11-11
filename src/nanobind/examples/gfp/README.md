# Approximate substructure searching with the formula fingerprint

The formula fingerprint counts atoms of each kind. It is very fast to compare, but it knows
nothing about connectivity, so it is not selective. Used with Tversky similarity it can find
the molecules that have at least the atoms of a query, so that substructure searching, or another
precise method, need only be applied to those. Its advantage is speed.

Run the scripts via the repository wrapper, from a virtual environment matching the installed
LillyMol Python ABI:

```shell
"$LILLYMOL_HOME/run_python.sh" "$LILLYMOL_HOME/src/nanobind/examples/gfp/formula_tversky.py"
```

| File | Purpose |
| --- | --- |
| `formula_tversky.py` | How the distances change with the Tversky weights, and how well the fingerprint works as a filter, with and without the counts of atoms bearing hydrogens, checked against real substructure searching. Uses `molecules.smi`, or `--file` and `--query`. |
| `formula_prefilter_timing.py` | Times the filter and the precise search on a large collection, checks that no molecule containing the query is rejected, and estimates the time saved. Needs a file made by `formula_fingerprint -H`. |
| `molecules.smi` | 40 molecules, with named examples of molecules that do and do not contain a phenol. |
| `molecules.formula_nh.gfp` | The fingerprints of those molecules, made with `formula_fingerprint -H molecules.smi`, for testing the timing script. Far too small for the times to mean anything. |

For a large collection:

```shell
formula_fingerprint -H collection.smi > collection.formula_nh.gfp
"$LILLYMOL_HOME/run_python.sh" formula_prefilter_timing.py collection.formula_nh.gfp --query 'Clc1ccccc1'
```

The hydrogen counts must be left out (`-H`, or `GFP.formula(hydrogens=False)`): with them,
molecules that contain the query are lost, because the query has hydrogens on atoms that are
substituted in the molecule.

## What to expect on a large collection

The measurements below are for the 2,748,104 molecules of ChEMBL 37, searched with
`Tversky(query=2, candidate=0)` at distance 0, on one thread. The precise search, which is
what the filter saves, is parsing the smiles and then a substructure search (about 6 plus 13
to 26 microseconds per molecule). Its times are estimated from a random sample of 20,000
molecules, each freshly parsed, since the first search of a molecule also perceives rings and
aromaticity.

| Query | Passes the filter | Substructure search of all | Filter, then search the rest | Times faster | Largest possible |
| --- | --- | --- | --- | --- | --- |
| benzene | 90.7% | 105 s | 96 s | 1.10 | 1.10 |
| phenol | 83.7% | 78 s | 66 s | 1.18 | 1.20 |
| pyridine | 60.9% | 79 s | 49 s | 1.63 | 1.64 |
| piperazine | 50.7% | 68 s | 35 s | 1.95 | 1.97 |
| chlorobenzene | 17.2% | 65 s | 11 s | 5.7 | 5.8 |
| trifluoromethyl | 8.8% | 69 s | 6 s | 11 | 11 |
| steroid core | 5.4% | 75 s | 4 s | 18 | 19 |

* The filter is not selective, so how much it saves depends on the query. For a generic query
  such as phenol it still saves a little.
* The saving can never be more than one over the fraction that passes, because the molecules that
  pass still have to be searched, however slow the precise method is. The measured speedups are
  within 1-2% of that limit, since the precise search costs hundreds of times more per molecule
  than the filter. A more expensive second stage, such as computing and comparing a richer
  fingerprint, cannot do better. Only a more selective query can.
* No molecule that contains the query was rejected, for any of the 11 queries, among the
  roughly 47,000 sampled molecules that contain one.
* Scanning the 2.7 million fingerprints takes about 75 ms, 27 ns per molecule, but asking for
  all the molecules that pass takes 0.1 to 0.6 s, because most of that time is building the list
  of results, so it grows with the number that pass.
* Loading the 2.7 million fingerprints, from a 345 MB file, takes about 2 seconds and
  0.6 GB of memory.

Scanning the same 500,000 molecules, with `Tversky(2, 0)` and one result wanted, costs about 25 ns
per molecule for the formula fingerprint without the hydrogen counts (33 ns with them),
48 ns for a stored path fingerprint (`GFP.iw()`) and 193 ns for a stored EC fingerprint
(`GFP.ec()`). The formula fingerprint is cheaper, but only by about a factor of two compared with
the path fingerprint, since there is a fixed cost per fingerprint. If a path fingerprint is
already stored for every molecule, filtering with the formula fingerprint first gains little. The
formula fingerprint is tiny (18 bytes), so it is most useful where the richer fingerprints, or a
substructure search, would have to be computed for each molecule that is examined.

Each script is also a Bazel test, named `//nanobind:gfp_formula_tversky_example_test`
and `//nanobind:gfp_formula_prefilter_timing_example_test`.
