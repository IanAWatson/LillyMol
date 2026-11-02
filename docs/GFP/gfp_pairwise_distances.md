# gfp_pairwise_distances

Compute distances for identifier pairs using a pool of LillyMol fingerprints:

```shell
gfp_pairwise_distances -p molecules.gfp pairs.txt > distances.txt
gfp_pairwise_distances -p molecules.gfp -K id1,id2 -o 3
```

Each input line starts with two identifiers. Output contains those identifiers
and their distance; `-a` also copies remaining input tokens. `-r N` instead
writes `N` random pairs with distinct identifiers.

Distances use `Fraction_as_String` for stable formatting across Linux and macOS.
The default is four decimal places, with trailing zeros removed. `-o N` changes
the number of decimal places, replacing the former significant-digit setting.
Higher precisions use the formatter without a lookup table to avoid excessive
memory use. Similarity calculations and distance filtering are unchanged.
