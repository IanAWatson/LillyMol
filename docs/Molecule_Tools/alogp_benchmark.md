# alogp_benchmark

Times LillyMol's ALogP (Wildman and Crippen) against RDKit's `Crippen.MolLogP` and
compares the values. It starts from a file of smiles and reports speed, correlation,
mean difference, the number of molecules that differ significantly, the most common
difference values with an example of each, and the largest outliers.

```
contrib/bin/alogp_benchmark.sh -o results molecules.smi
```

The input is `smiles identifier` per line. `-n 10000` uses only the first 10,000
molecules, `-t` sets the difference counted as significant (default 0.1), and `-k` how
many outliers to list. It writes `results/values.tsv` (one column per method) and
`results/results.json`.

It needs `rdkit` and `numpy` in addition to the LillyMol python bindings, so the
wrapper runs through `run_python.sh`; set `PYTHON` to an interpreter that matches
`lib/lillymol.so.soabi` and has rdkit installed.

## What is measured

| | |
|---|---|
| RDKit | `MolFromSmiles`, then `MolLogP` (which includes adding hydrogens) |
| LillyMol python | `Molecule.build_from_smiles`, then `ALogP.logp` |
| LillyMol command line | one run of `alogp` over the whole file, parsing included |

Each LillyMol method is run twice, with default settings and with the RDKit
compatibility settings (`-Y alcacid -Y RDKIT.N+ -Y RDKIT.HP` on the command line,
`set_use_alcohol_for_acid` and `set_rdkit_phoshoric_acid_hydrogen` in python). Timings
are single threaded wall clock from one run, so treat the ratios as approximate.

The default settings follow the paper, so carboxylic acids (O12) differ from RDKit,
which types them as alcohols; with 100,000 ChEMBL molecules this is the main source of
the default-setting difference. Use the compatibility settings to look for the real
differences.

## Testing a new build

By default the command line tool is `$LILLYMOL_HOME/bin/Linux/alogp` and the bindings
are those in `$LILLYMOL_HOME/lib`. To benchmark a fresh build without replacing the
installed one, stage the bindings in a scratch directory and point at the new tool:

```
src/copy_shared_libraries.sh /tmp/newlib
PYTHONPATH=/tmp/newlib LD_LIBRARY_PATH=/tmp/newlib \
  ALOGP=src/bazel-bin/Molecule_Tools/alogp \
  python3 contrib/bin/alogp_benchmark.py -o results molecules.smi
```

## Example

On 100,000 ChEMBL molecules (single thread):

| | seconds |
|---|---|
| RDKit parse + MolLogP | 10.2 + 24.0 |
| LillyMol python parse + logp | 0.41 + 1.69 |
| LillyMol command line, whole run | 1.83 |

The logp step is about 14 times faster, and parse plus logp 16 to 18 times faster. With
the compatibility settings the Pearson correlation with RDKit is 0.9997, the mean
absolute difference 0.002, and 228 molecules differ by more than 0.1. Those remaining
differences are mostly places where LillyMol's aromaticity perception differs from
RDKit's, and quaternary aryl ammoniums, where RDKit's N13 pattern does not match and it
assigns the catch-all N type.
