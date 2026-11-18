# mkaspirin

`mkaspirin.py` generates deliberately complicated recursive SMARTS for testing
query parsers and matchers. It starts from a fixed aspirin query and repeatedly
replaces an atom description with redundant rooted recursive context. Every
atom token retains its source aspirin atom index, so a carbonyl oxygen is
expanded using context rooted at that oxygen, rather than at the methyl carbon
where the usual aspirin SMILES starts.

The existing `contrib/bin/mkbenzene.py` remains available for the simpler,
symmetric benzene case. This is a separate generator; it uses standard-library
`argparse` and `random`, with LillyMol bindings for molecular perception.

## Generate a large test query

Run through the wrapper with a Python interpreter matching the staged LillyMol
bindings (`lib/lillymol.so.soabi`):

```shell
PYTHON=/path/to/venv/bin/python3 \
  ${LILLYMOL_HOME}/contrib/bin/mkaspirin.sh \
  --target_bytes 3000000 --write_last --seed 19 > aspirin.smt
```

The wrapper calls `run_python.sh`. Python 3.10 or newer and staged nanobind
bindings are required. No NumPy, absl or RDKit is needed to generate queries.

Each output line contains one SMARTS and a name such as `#137`. The name is the
number of expansions, not part of the SMARTS. Diagnostics go to stderr. Test it
with LillyMol:

```shell
printf '%s\n' 'CC(=O)Oc1ccccc1C(=O)O aspirin' > aspirin.smi
tsubstructure -q S:aspirin.smt -m matches aspirin.smi
```

Portable SMARTS is the default. To include LillyMol extensions:

```shell
${LILLYMOL_HOME}/contrib/bin/mkaspirin.sh \
  --target_bytes 3000000 --write_last --seed 19 \
  --lillymol_extensions > aspirin_lillymol.smt
```

Ring-atom predicates then include `/IWfss1`, `/IWfss>0` or `/IWfss<2`, based on
aspirin's one-ring fused system. Relational atomic predicates such as `D>1`,
`D<3`, `H<2`, `r>5` and `r<7` can also be selected. The script places `/IW`
directives at the front of atomic predicates by convention. These queries are
for LillyMol; they are not intended for RDKit or another standard-SMARTS parser.
`--no_lillymol_extensions` explicitly selects portable mode.

## What remains invariant

The baseline is a strict query derived from `CC(=O)Oc1ccccc1C(=O)O`. Every atom
has constraints for element, aromaticity, hydrogen count, degree, total
connectivity, ring membership, ring-bond count, isotope zero and charge zero.
Ring atoms also require a six-membered ring. Bonds are explicit: aromaticity is
checked before the underlying Kekule single/double type.

Every generated atom retains those baseline constraints. Additional atomic
predicates repeat facts already implied by them. Recursive contexts are source
subgraphs rooted at the same source atom, so any baseline embedding also
supplies an embedding of the added context. The generator therefore aims to
preserve the baseline match set, rather than merely guaranteeing that aspirin
itself matches. It does not generate weakened or intentionally mismatching
queries.

This is a substructure query, not a molecular-identity function. Disconnected
additional components do not prevent a match. Aspirin anion, an isotope-labelled
aspirin, positional isomers, methyl ester and other related structures can fail
the strict baseline. Their failure is not evidence of a matcher bug.

## Controls

| Option | Meaning |
| --- | --- |
| `--nsmarts N` | Perform N expansions; default 100. With zero, write the baseline once. |
| `--target_bytes N` | Expand until the SMARTS reaches at least N ASCII bytes, excluding its name and newline. Cannot be combined with `--nsmarts`. |
| `--write_last` | Write only the final expression. Recommended for large runs. |
| `--seed N` | Reproduce the random choices; default 2026. |
| `--max_depth N` | Maximum number of nested `$()` levels; default 32. |
| `--alternatives N` | Add N ORed recursive contexts at each expansion; default 1. |
| `--contexts LIST` | Choose from `atom,bond,radius2,molecule`; default `bond,radius2,molecule`. |
| `--lillymol_extensions` | Add fused-system and relational LillyMol constraints. |
| `--no_lillymol_extensions` | Explicitly select portable SMARTS, the default. |

`atom` repeats only the root atom. `bond` includes all atoms within one bond of
the root; `radius2` includes those within two bonds. `molecule` includes the
whole aspirin graph. Each new context gets independently chosen atom predicates
and ring labels. The first atom in every context is the atom being expanded.

For example, generate OR alternatives from full-molecule contexts:

```shell
${LILLYMOL_HOME}/contrib/bin/mkaspirin.sh \
  --nsmarts 100 --contexts molecule --alternatives 2 --write_last > aspirin.smt
```

A size target can overshoot by one expansion. Depth is a hard limit: if no
expandable atoms remain before the requested count or size is reached, the
script reports an error and exits nonzero. Without `--write_last`, previously
written expressions remain on stdout. Increasing `--max_depth` controls the
generator; it does not guarantee that a matcher can safely parse or search that
depth. Matcher stack, memory and runtime limits still apply.

## Implementation and validation

The source graph is perceived once with LillyMol. Small rooted query templates
are generated from that graph while retaining atom identities. Expanded
expressions are stored as structured nodes. Random expansion sites are kept in
a swap/pop list, avoiding repeated insertion into a growing flat token array.
The SMARTS size is maintained incrementally, and serialization is iterative.

With `--write_last`, the growing expression is serialized once at the end.
Without it, serializing every intermediate expression can dominate runtime and
produce a very large output file. The expansion sequence for a given seed is
the same in both cases. Reproducibility is for the same script/runtime and
LillyMol perception; it is not a promise across future algorithm changes.

Run the focused tests using staged bindings:

```shell
PYTHON=/path/to/venv/bin/python3 \
  ${LILLYMOL_HOME}/run_python.sh \
  ${LILLYMOL_HOME}/contrib/bin/test/mkaspirin_test.py
```

They check all 13 root atoms at all four context sizes, multiple seeds and OR
alternatives, portable and extension modes, randomized aspirin SMILES, a panel
of related structures, byte accounting, depth exhaustion, CLI validation and
serialization deeper than Python's usual recursion limit. When RDKit is
installed, it independently parses and searches representative portable queries.
No test requires RDKit to understand LillyMol extensions.

Portable and extension-mode expressions of at least 3 MB have also been generated
and matched against aspirin with `tsubstructure`. This is a stress-test scale
check, not a general performance guarantee for arbitrary molecules or matchers.


## Python-binding benchmark

[Recursive aspirin SMARTS: Python binding timings](mkaspirin_benchmark.md)
compares LillyMol and RDKit parsing, first and repeated searches, mixed-panel
searches and peak process memory from 1 KB to 3 MB. It includes raw measurements,
hardware and methodology, and a driver for repeating the experiment.
