# gfp_spread_buckets

Selects diverse molecules while balancing selections across buckets. Each
fingerprint must contain integer bucket values in the specified tags.

```shell
gfp_spread_buckets -B NATOMS -n 100 molecules.gfp > selected.tdt
```

Repeat `-B` to balance more than one bucket classification. `-o` requests strict
rotation across buckets. `-A previous.gfp` supplies previously selected
molecules; add `-a` when that file also contains bucket information. `-p WEIGHT`
reads a distance scale factor from the specified tag.

Output `DIST` values use platform-independent formatting with four decimal
places and trailing zeros omitted, including scaled distances above one.
Selection calculations retain their original precision. The first selection
has no `DIST` when it has no nearest selected neighbour.

