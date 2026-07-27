# chemotypes

## Overview

`chemotypes` extracts query-defined ring-system chemotypes from molecules. A
substructure query identifies the seed ring system, and the output molecule is
reduced to that ring system plus any optional neighbouring ring systems,
linkers, exocyclic double-bonded atoms, or attachment-point annotations.

The intended use is grouping molecules by a chemically meaningful ring-defined
core. Chemists often describe a series by a chemotype such as a piperidine,
benzimidazole, or substituted fused ring system. This tool makes that
definition explicit by combining a query with deterministic rules for expanding
from the matched seed ring system.

## Quick Start

Extract the ring system containing the first atom matched by a SMARTS:

```shell
chemotypes -s '[N]' -z i -i smi input.smi > chemotypes.smi
```

For example, from:

```text
CN1CCCCC1 mol1
CC miss
```

the command writes:

```text
N1CCCCC1 mol1
```

By default, a molecule that does not match any query is an error. Use `-z i`
or `-z ignore` to ignore nonmatching molecules.

## Query Matching

The seed must be specified with one or more queries:

```shell
chemotypes -q seed.qry input.smi > chemotypes.smi
chemotypes -s '[nH]1cccc1' input.smi > chemotypes.smi
```

Queries are examined in command-line order. The first matching query is used.
By default, that query must identify exactly one ring system across all of its
embeddings. This means a broad query such as `c` is acceptable for a molecule
with one aromatic ring system, but ambiguous for biphenyl because the embeddings
touch two different ring systems.

Within the accepted query, the first embedding that touches the identified ring
system is used. The first matched atom in that embedding that belongs to the
ring system defines the seed atom.

Use `-z f` or `-z first` to keep the historical behaviour of using the first
embedding of the first matching query, even if the query touches multiple ring
systems.

The query should normally include ring constraints where appropriate. If the
first matching query has no matched atom in a ring system, processing stops for
that molecule.

## Adjacent Ring Systems

Use `-n <n>` to include directly adjacent ring systems.

```shell
chemotypes -s '[c]' -n 1 -i smi input.smi > chemotypes.smi
```

Ring systems are considered directly adjacent when they can be reached from the
seed ring system through a path whose internal atoms are not ring atoms. Ring
systems beyond another intervening ring system are not considered directly
adjacent.

The adjacent systems are ranked by distance from the seed ring system, with
query atom order used to make deterministic choices in ambiguous cases. If a
molecule has more adjacent systems than requested, only the best ranked systems
are included.

Use `-t` with `-n` to include all systems tied at the cutoff distance:

```shell
chemotypes -s '[c]' -n 1 -t -i smi input.smi > chemotypes.smi
```

In a case where two phenyl rings are directly attached to the seed ring at the
same distance, `-n 1` includes one of them, while `-n 1 -t` includes both.

Use `-D <dist>` to limit adjacent ring systems by maximum bond distance from
the seed ring system.

## Linkers and Exocyclic Atoms

When adjacent ring systems are included, linker atoms between selected ring
systems are included by default. The linker is the non-ring connected component
that connects the selected ring systems, not just one shortest path.

Terminal exocyclic atoms double bonded to retained atoms are also included by
default. This preserves common ring-system features such as carbonyl oxygen
atoms.

## Attachment Points

The default output contains only the retained chemotype core and does not encode
how that core was attached to the rest of the parent molecule.

Use `-I <iso>` to label retained ring atoms that have an external attachment:

```shell
chemotypes -s '[N]' -I 99 -i smi input.smi > labelled.smi
```

This labels the ring atom exit point with the fixed isotope value. `-I` cannot
be used with atom typing.

Use `-u` to retain one-hop atoms attached to retained ring atoms without
applying atom-type labels:

```shell
chemotypes -s '[N]' -u -i smi input.smi > attached.smi
```

When `-P` is specified, atom typing is applied before the molecule is reduced.
The tool automatically retains one-hop atoms attached to the chemotype core and
keeps isotopic atom-type labels on non-terminal attachment atoms. Isotopes on
the core chemotype atoms and on atoms that were singly connected in the parent
molecule are cleared. This can distinguish, for example, the same ring system
attached through carbon or nitrogen.

```shell
chemotypes -s '[N]' -P UST:ARY -i smi input.smi > typed.smi
```

Use `-x` with `-u` to ignore singly connected attached atoms when retaining
untyped attachment atoms:

```shell
chemotypes -s '[N]' -u -x -i smi input.smi > attached.smi
```

## Parent Output

Use `-p <text>` to write the parent molecule before each generated chemotype:

```shell
chemotypes -s '[N]' -p PARENT -z i -i smi input.smi > parent_and_chemotype.smi
```

The value is appended to the parent molecule name as a separate token, which
makes parent records easier to distinguish from chemotype records when reviewing
SMILES output. Use `-p .`, `-p def`, or `-p default` to write the parent without
adding an annotation.

This is mainly useful for debugging, review, and explaining why a chemotype was
generated. The parent is written only for molecules that successfully generate a
chemotype.

## Summary Output

Use `-F <fname>` to accumulate generated chemotypes as
`dicer_data::DicerFragment` textprotos:

```shell
chemotypes -s '[N]' -z i -i smi -F summary.textproto input.smi > chemotypes.smi
```

The summary is keyed by the unique SMILES of each generated chemotype. The first
occurrence creates a proto with:

| Field | Meaning |
| --- | --- |
| `smi` | Unique SMILES of the chemotype. |
| `par` | Name of the first parent molecule that generated the chemotype. |
| `nat` | Number of atoms in the chemotype. |
| `n` | Number of parent molecules that generated this chemotype. |
| `iso` | Present when isotopic labels encode attachment information. |

For example, if two input molecules reduce to the same piperidine chemotype,
the summary contains:

```text
smi: "N1CCCCC1" par: "mol1" nat: 6 n: 2
```

When `-P` is active, `iso` is written as `ATYPE`. When `-I` is active, `iso` is
written as `ATT`.

## Filtering by Ring Count

Use `-r <n>` to process only molecules with at least `n` rings:

```shell
chemotypes -s '[c]' -n 1 -r 2 -z i -i smi input.smi > chemotypes.smi
```

This is useful when the requested chemotype expansion requires nearby rings and
single-ring molecules should simply be skipped.

## Options

| Option | Meaning |
| --- | --- |
| `-q <query>` | LillyMol query file defining the chemotype seed. May be repeated. |
| `-s <smarts>` | SMARTS defining the chemotype seed. May be repeated. |
| `-n <n>` | Include `n` directly adjacent ring systems. Default `0`. |
| `-t` | With `-n`, include all adjacent ring systems tied at the cutoff distance. |
| `-D <dist>` | Maximum bond distance to an adjacent ring system. |
| `-r <n>` | Minimum number of rings required for processing. |
| `-u` | Include one-hop atoms attached to retained ring atoms. |
| `-x` | With `-u`, ignore singly connected attached atoms. |
| `-I <iso>` | Label retained ring exit-point atoms with isotope `iso`. Incompatible with `-P`. |
| `-P <atype>` | Atom typing specification for retained non-terminal attachment atoms. |
| `-p <text>` | Write the parent molecule before each generated chemotype; append `<text>` unless it is `.`, `def`, or `default`. |
| `-F <fname>` | Write accumulated `dicer_data::DicerFragment` textproto summary. |
| `-z i`, `-z ignore` | Ignore molecules that do not match any query. |
| `-z f`, `-z first` | Use the first embedding when a query matches multiple ring systems. |
| `-S <stem>` | Output file stem. Default stdout. |
| `-i <type>` | Input type. Usually inferred from filename. |
| `-o <type>` | Output type. Default SMILES. |
| `-g ...` | Chemical standardisation. |
| `-l` | Reduce to largest fragment. |
| `-c` | Remove chirality. |
| `-A ...` | Aromaticity options. |
| `-E ...` | Element options. |
| `-v` | Verbose progress and summary reporting. |

