# gfp_profile_activity_by_bits

Profile activities for molecules with and without each fingerprint bit:

```shell
gfp_profile_activity_by_bits -E activity.txt molecules.gfp > profiles.txt
gfp_profile_activity_by_bits -E activity.txt -M prec=6 -g -j molecules.gfp
```

Result-file floating-point values use `Fraction_as_String` for stable formatting
across platforms. The default is four decimal places, with trailing zeros
removed; `-M prec=N` changes the number of decimal places. This replaces the
former significant-digit formatting. Signed values and values outside `[0,1]`
are supported without allocating a lookup table.

The formatter covers fixed/sparse fingerprint statistics, classification `NB`
log-ratios (`-C`), and statistical table output (`-T filename`). Integer counts,
bit identifiers, missing-value markers, and calculations are unchanged.

Classification runs now return zero on success; earlier versions incorrectly
returned one even after producing valid output. Read and output failures still
return a nonzero status.
