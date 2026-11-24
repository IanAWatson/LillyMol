# iwstats

`iwstats` is a command line tool that computes a variety of statistical
measures comparing an experimental set of continuous values and a
predicted set of such values.

The observed and predicted values can both be in the same file, or can come
from different files.

Note that this could obviously be done using R or Julia or Python but in
an environment where pre-requisites have often been difficult to achieve
this works well and is fast. And some measures have arisen at Lilly, and
this is the reference implementation.

## Minimal Invocation.
If you have a tabular file with the experimental values in column 2
and the predicted values in column 3
```
iwstats -e 2 -p 3 file.dat
```
will yield various statistical association measures between those two columns. If
the input file has a header record, add the `-j` option, or `-s 1`.

If the experimental results are in a different tabular file, `id expt` form
```
iwstats -E expt.dat -e 2 -p 2 -j file.dat
```

## Options
The tool is complex and recognises the following options
```
Computes Bsquared and other statistics - allows missing values
 -e <col>       column for experimental (measured) values
 -E <file>      activities are in a different file <file>
 -z             strip leading 0's from identifiers when using -E
 -p <col>       column for predicted values
 -s <number>    skip <number> records at the top of each file
 -n <number>    process only the first/best <number> records
 -t <number>    compute Bsquared for sections of the data. For example,
                -t 10 would report values for the first/best 10, 20, 30 values
 -P <number>    sample the first/best <number> percent of the data
 -j             treat as a descriptor file
 -M <string>    missing value string
 -q             quietly ignore missing values
 -R .           randomise the sort when duplicate predicted values present
 -c <number>    number of valid pairs needed for producing correlations (default 20)
 -w             when duplicate values present, suppress computation of best
                and worst BSquared values
 -r <float>     max relative error allowed computations
 -T             truncate predicted values to the experimental range
 -h             use traditional Q2 formula
 -b <nbucket>   compute distribution functions across <nbucket> buckets
 -F <f>         calculate number of predictions outside <f> fold of experimental (multiplicative)
 -D <f>         calculate number of predictions differing by <d> from experimental (additive)
 -k             just skip predicted values that have no experimental value (repeat for quiet)
 -d             compute Dave's cMSD rather than cMSR
 -o <float>     cutoff for active/inactive (for enrichment metrics BEDROC EF ...)
 -a <float>     BEDROC alpha value (default 20)
 -f <float>     Enrichment Factor default fraction (default 0.5)
 -u <float>     discard experimental values below <float> - useful for studying actives only 
 -U <fraction>  keep only the <fraction> most active experimental values
 -m <col/name>  do analysis by data in column <col> or descriptor name <name>
 -L <fname>     write residuals to <fname>
 -v             verbose output
```

Some combination of `-E`, `-e` and `-p` must be used in order to specify where the
column of experimental data, and the column of predicted values are located.

### Missing Values
Within LillyMol many tools use `.` as the missing value string. That can be changed
with the `-M` option. When missing values are present, various warnings will be
issued, but if the `-q` option is used, most of those warnings will not be issued.

## Precision and valid pairs

Measured and predicted values, residuals, and working arrays use double precision.
Correlations and correction coefficients use centered sums to avoid subtracting
nearly equal raw sums for values with large offsets.

The default missing-value token is `.`; `-M` selects another token. Rows missing a
measured value are excluded. A missing prediction excludes that pair only from
the corresponding prediction column. Measured-value filters (`-u` and `-U`) are
applied before ranking and subsetting. Counts, errors, correlation, enrichment,
and fold-error fractions use retained valid pairs. With `-T`, truncation affects
predicted summaries and error calculations; the raw prediction ranking is
retained. `-B 100` reports the maximum absolute error. `-W` works independently
of `-w`.

`R2` is the **signed square of Pearson correlation**, rather than the regression
coefficient of determination. For conventional predictive R-squared, use `-h`
and read `trQ2`, which is `1 - SSE/SST`. Without `-h`, `Q2` uses the historical
bounded transformation when SSE exceeds SST. Constant vectors report correlation
zero while retaining their error statistics. When measured values are constant,
Q2 is 1 for perfect predictions and 0 otherwise. Zero residual spread is handled
without dividing by zero in cMSR/cMSD calculations.

By default, predictions tied at the same value retain input order for the input
Bsquared calculation. `-R .` explicitly randomizes ties. Best/worst Bsquared still
bounds the possible orderings of ties unless `-w` suppresses that calculation.
Replacing the previous sorter can therefore change `B2.input` on tied data.

`-H` coverage and `-C` correction outputs require one input file, one predicted
column, and no subsets or marker groups. Ambiguous combinations are rejected
rather than overwriting earlier results. Coverage requires more than 100 valid
pairs. Read errors and detected output write/flush failures produce nonzero exit
status. Residual output is explicitly flushed and closed at completion.
