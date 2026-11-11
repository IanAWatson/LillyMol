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

How much the filter saves depends on the query, because it is not selective. On 200,000
ChEMBL molecules, estimated from a sample of 20,000, a phenol query passes 84% of the molecules
and the filter saves little (1.2 times faster than searching them all), whereas chlorobenzene
passes 15% (6.5 times faster) and a sulfonamide passes 18% (5.5 times faster).
The hydrogen counts must be left out (`-H`, or `GFP.formula(hydrogens=False)`): with them,
molecules that contain the query are lost, because the query has hydrogens on atoms that are
substituted in the molecule.

Each script is also a Bazel test, named `//nanobind:gfp_formula_tversky_example_test`
and `//nanobind:gfp_formula_prefilter_timing_example_test`.
