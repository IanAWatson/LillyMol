# substitutions

`substitutions` counts bonds from query-matched atoms to unmatched atoms and can
extract the substituent fragments. It can also write an R-group table containing
fragment SMILES and optional molecular properties.

```shell
substitutions -s 'N1CCCCC1' -z first input.smi > counts.txt
substitutions -s 'N1CCCCC1' -z first -x -S fragments input.smi
substitutions -s 'N1CCCCC1' -z first -x -w -d groups.tsv input.smi
```

The first command writes `Name SUBS0` followed by a substitution count for each
molecule. `-S fragments` writes `fragments.smi`; `-x` excludes the matched anchor
atom from each fragment. Parents are included in `-S` output by default; `-y`
suppresses them. `-d` names a tab-separated R-group file. `-w` adds properties to
that table; it does not change the normal count output. `-n` suppresses counts.

## Queries and attachment points

Use repeated `-s` SMARTS or `-q` query files for multiple query definitions.
`-e` builds queries from SMILES. `-m` instead reads a single file of isotopically
labelled query molecules; its isotopes define attachment points. This form must
not be combined with `-s`, `-q`, `-e`, or another `-m` file. `-h` allows more than
two isotopic attachment points in `-m` input.

For fragment extraction, `-j` specifies zero-based matched query atom indices:

```shell
substitutions -s 'N-C-O' -j 0,2 -x -S fragments input.smi
substitutions -s 'N.N' -j 0-1 -z first -x -S linkers input.smi
```

A comma list specifies independent sites. A pair such as `0-1` requests the
fragment between the two sites; a one-atom linker is retained. Repeated `-j`
values correspond to successive queries. Without `-j`, all matched query atoms
are considered. The normal `SUBS` columns count outgoing bonds over the complete
retained embedding, rather than only the `-j` sites. An invalid attachment index
causes failure instead of an out-of-bounds access.

## Multiple matches

- `-z each` processes all embeddings.
- `-z first` retains the embedding with the most outgoing substitution bonds,
  using the selected attachment points where specified; ties retain the first.
  Despite its historical name, this is not always the first initial embedding.
- `-z ignmmatch` ignores queries with multiple retained matches.
- `-z ignore` allows individual queries not to match.
- `-z oknomatch` also allows a molecule to match none of the queries.

Without an explicit multiple-match policy, unresolved multiple matches cause
failure. Duplicate elimination compares matched atom sets when no attachment
points are specified, or corresponding selected sites for single-site lists.
Distinct matches are not silently treated as duplicates. Paired attachment
specifications bypass this duplicate-elimination step.

## R-group tables and missing fragments

`-d` uses one internal slot per matched query atom. R-group columns are retained
for positions substituted in at least one row for the same query. Unsubstituted
positions are represented by `H` (or the configured capping group). Missing or
unavailable fragments, including fragments rejected by `-C`, are represented by
`.` with missing property values. Such entries do not shift later groups.
Rows from queries exposing fewer groups are padded to the common header width.
A query with no substituted positions can produce a scaffold-only row.

Fragments spanning multiple attachment positions are not assigned to a single
R-group column: their affected slots remain `.`. Use fragment output (`-S` or
`-J`) when those fragments are needed. `-r` permits ring substituents, and `-c`
labels atoms where the traversal reconnects to the matched query. Temporary
attachment and ring-closure labels are restored on the parent before later
matches and parent output.

`-C n` limits fragment size during traversal. The single-site traversal counts
the anchor before `-x` removes it; pair-based traversal counts interior atoms.
`-P file` supplies a capping group and changes the unsubstituted-group reference.

## Properties and TPSA

The R-group properties are `natoms`, `nrings`, `amw`, `ro5_ohnh`, `ro5_on`,
`nvrtspsa`, `halogen`, `fcsp3`, `ringsys`, `arring`, `alring`, `ringatom`, and
`mxdst`. Molecular weight ignores isotope labels. `fcsp3` is the number of
saturated carbon atoms divided by the total atom count. `mxdst` describes the
largest connected fragment; calculating it does not remove components from the
stored group. Properties describe the extracted group, including retained
anchors or caps when requested.

`nvrtspsa` now uses `nvrtspsa::NovartisPolarSurfaceArea`, shared with the other
LillyMol descriptor tools. This intentionally replaces the old private table:

- aromatic sulfur and thioether sulfur receive the library's nonzero
  contributions (28.24 and 25.30 for the usual neutral types);
- explicit hydrogens are handled on a hydrogen-suppressed copy for TPSA;
- unclassified atom types follow the library's fallback contribution policy
  rather than always contributing zero.

Consequently, some TPSA values change. This is the default LillyMol calculation,
not RDKit compatibility mode. Other property definitions are retained.

## Output safety and errors

`-S`, `-Z`, `-J`, `-T`, and `-d` must not overwrite molecular input files. The
text-output checks include symbolic and hard-link aliases. `-T` appends `.smi`
when needed before checking. Unreadable input, invalid requested attachment
indices, and an unavailable R-group output file produce nonzero exit status.
If input processing fails, a partial R-group table is not emitted.

## Tests

From `src`:

```shell
bazel test -c opt --config=native //Molecule_Tools:substitutions_test
ruby test/run_all_test.rb -bindir "$PWD/bazel-bin/Molecule_Tools" substitutions
```
