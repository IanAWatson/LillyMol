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
