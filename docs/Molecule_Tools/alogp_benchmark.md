# alogp_benchmark

Times LillyMol's ALogP (Wildman and Crippen) against RDKit's `Crippen.MolLogP` and
compares the values. It starts from a file of smiles and reports speed, correlation,
mean difference, the number of molecules that differ significantly, the most common
difference values with an example of each, and the largest outliers.

```
contrib/bin/alogp_benchmark.sh -o results molecules.smi
```

The input is `smiles identifier` per line. `-n 10000` uses only the first 10,000
molecules and `-s 5000` skips the first 5,000 first (the tool is given the same range with
`-i do=` and `-i skip=`, so no temporary file is written). `-t` sets the difference counted as significant (default 0.1), and `-k` how
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

Output for 100,000 ChEMBL molecules, single threaded (RDKit 2026.09.1; RDKit failed to
parse 2 of them, which are left out of the comparison):

```
Timings (seconds, single thread)
  RDKit        parse    10.20  logp    24.05  (2 failed to parse)
  LillyMol py  parse     0.59  logp     1.69  (rdkit compat 1.95)
  LillyMol CLI whole run 1.83  (rdkit compat 1.83)
  speed up: logp only 14.2x, parse+logp python 15.0x, CLI 18.7x

Agreement with RDKit (LillyMol - RDKit); threshold 0.1
                           n  pearson     MAE     bias    RMSE   >thr     >1     max
  py_default           99998  0.98074  0.1098  -0.1091  0.3742   9693   9046   6.123
  py_rdkit_compat      99998  0.99971  0.0017  -0.0010  0.0426    228     82   3.600
  cli_default          99998  0.98074  0.1100  -0.1091  0.3742   9693   9046   6.123
  cli_rdkit_compat     99998  0.99971  0.0019  -0.0010  0.0426    228     82   3.600
```

The logp step is about 14 times faster than RDKit, and parse plus logp 15 to 19 times
faster. The Python and command line values agree to the three decimals the tool
prints. The python parse time varies from run to run (0.41 to 0.59 s in repeated runs);
the logp time is stable to a few percent.

With the compatibility settings, 228 molecules (0.2%) differ from RDKit by more than
0.1. Most of what remains is not an ALogP issue:

| difference | molecules | cause |
|---|---|---|
| -1.05, -0.97, -0.69, -0.60, -0.33 and similar | about 90 | LillyMol perceives more atoms as aromatic in fused ring systems that contain a ring N-C(=O)-N, so different atom types are assigned |
| +0.14 | 22 | aryl quaternary ammonium, `Ar-N+(C)(C)C`: RDKit's N13 pattern needs aliphatic neighbours, so it uses the catch-all N type; LillyMol uses N13 |
| +0.05 | about 200 (below the 0.1 threshold) | quinoidal tautomer drawings, where RDKit leaves the ring fusion bond non-aromatic |
| the rest | | malformed or charge-separated drawings, and a few cationic dyes |

In default mode the large count (9,693) is almost entirely carboxylic acids, where the
paper assigns O12 and RDKit assigns O2.
