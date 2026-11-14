# gfp_lnearneighbours

Find neighbours by loading target fingerprints into memory and scanning a
fingerprint haystack:

```shell
gfp_lnearneighbours -p targets.gfp -n 5 -h haystack.gfp > neighbours.nn
```

Text/TDT distances use `Fraction_as_String` with three decimal places and
trailing zeros removed, avoiding platform-dependent scientific notation.
This replaces the former three-significant-digit formatting. Neighbour searches,
ordering, thresholds, and counts are unchanged.

The writer honours `-D tag` and `-B DIST=tag` for distance tags and `-B nosmiles`
for omitting SMILES. `-k` writes NearNeighbours textproto records instead;
protobuf distances retain their original numeric precision.

## Sorting by proximity to reference molecules (DISE-style ordering)

For Directed Sphere Exclusion (DISE)-style prioritisation, first measure how
close each candidate is to a reference set. Put the candidates in the `-p`
file so the output contains one record per candidate:

```shell
gfp_lnearneighbours -p candidates.gfp -n 1 references.gfp > file.nn
tdt_sort -T DIST file.nn > file.sorted.nn
```

The default nearest-neighbour output is TDT (Daylight Thor DataTree), an
obsolete but useful record-based file format. The nearest neighbour's distance
is stored in the first `DIST` tag. `-n 1` retains just that neighbour;
`tdt_sort -T DIST` sorts numerically on the first occurrence, placing the
closest candidates first. Whole TDT records are retained together.

To reverse the ranking and put candidates farthest from the reference set first:

```shell
tdt_sort -T DIST -r file.nn > file.sorted.nn
```

This ranks each candidate by its minimum distance to any reference molecule.
The original DISE reference ranking instead uses similarity to the first
reference, then subsequent references to break ties. Both provide input
priorities for [leader/sphere-exclusion selection](gfp_leader.md#dise-style-input-ordering).
Use the candidate ranking to order molecules before generating their leader
input fingerprints. Sorting the neighbour report preserves neighbour records,
including their identifiers and distances.
