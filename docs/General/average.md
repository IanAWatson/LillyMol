# average

`average` reports counts, minimum and maximum values, sums, averages and sample
standard deviations for selected numeric columns.

```shell
average -c 2 data.txt
average -d activity descriptors.txt
```

Column numbers are one-based. `-d` selects columns from a header; `-j` skips a
header when columns are selected numerically. Multiple input files are reported
independently.

## Prevalence-weighted statistics

Use `-h` to specify a column containing the number of observations represented by
each row. For example, with `counts.txt` containing:

```text
1 1
9 3
```

```shell
average -c 1 -h 2 -p 9 -e 0,50,100 counts.txt
```

This represents four observations: `1, 9, 9, 9`. The average is 7; there are three
instances of 9, with fraction 0.75. The reported percentiles are 1, 9 and 9.
Zero-prevalence rows contribute nothing, including to percentiles and special
value counts. An input containing only zero-prevalence rows produces no
statistics and exits successfully.

Prevalences must be non-negative integers. `-h mult=10`, in addition to a column
specification, scales fractional prevalences to integer counts. Overflow is
rejected.

A prevalence column can also be specified by name:

```shell
average -d activity -h count first.txt second.txt
```

Named columns require a header (`-d`, `-j` or `-s`). Names are resolved separately
for each file, so columns can move between files. Missing or duplicate prevalence
names are errors. Numeric column specifications and multipliers apply to every
input file.

## Percentiles and memory

`-e` accepts individual percentiles or comma-separated lists, including 0 and
100. After sorting the represented observations, the tool selects the nearest
index to `percentile * (number_of_observations - 1) / 100`; it does not interpolate
between observations.

Without `-r`, percentiles retain all input values. With prevalence weighting,
only one value/count pair is stored per valid input row, regardless of the count.
Weighted percentiles therefore describe the full represented population without
expanding repeated observations in memory.

For unweighted input, `-r N` bounds storage to a reservoir sample of at most N
values (N must be at least 10). Only percentiles use the sample; the other
statistics use all observations. Sampled percentiles are approximate and can vary
between runs. Combining `-r` with `-h` is rejected; omit `-r` for exact weighted
percentiles.
