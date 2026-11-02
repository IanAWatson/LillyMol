# gfp_distance_matrix

Computes a distance matrix from LillyMol fingerprints.

```shell
gfp_distance_matrix -h molecules.gfp > distances.txt
gfp_distance_matrix -m molecules.gfp > similarities.txt
gfp_distance_matrix -l molecules.gfp > lower_triangle.txt
```

Normal text output uses platform-independent formatting with four decimal
places and trailing zeros omitted. Distances are calculated at their original
precision. `-g tab` selects a tab separator; `-h` adds a header.

`-C other.gfp` compares fingerprints against another set. `-e` reduces memory
use by recomputing distances rather than storing the symmetric matrix.

The existing `-f` mode uses aligned two-decimal strings; repeat `-f` to squeeze
multiple blanks. This mode retains its original formatting and ignores `-g`.
