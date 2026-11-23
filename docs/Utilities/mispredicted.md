# mispredicted

`mispredicted` combines prediction records from multiple files and identifies
molecules with large average prediction errors. These can be candidates for
investigating activity cliffs or measurement problems. The tool reports prediction
performance; it does not determine why a molecule is poorly predicted.

```shell
mispredicted -E measured.txt -n 2 -w 100 -d pred1.txt pred2.txt pred3.txt
```

This writes at most 100 identifiers with at least two valid prediction records,
ranked by mean absolute error. Identifiers below the minimum are excluded before
ranking and selection.

## Input files

Experimental and prediction files require a header. The first column is the
identifier; the second is the measurement or prediction by default:

```text
ID measured
mol1 2.0
mol2 4.0
```

```text
ID predicted
mol1 3.0
mol2 3.5
```

Use `-e <column>` and `-p <column>` to select other experimental and prediction
columns. Column numbers are one-based and must be at least 2. Continuous input
should use spaces as separators; classification input also translates tabs.
Identifiers are matched exactly and must be single tokens. Measurements and
predictions must use the same activity scale.

Repeated experimental identifiers are errors unless `-h` is specified; `-h`
retains the first value. Prediction files may contain different subsets of the
experimental identifiers. `-g` ignores predictions with no experimental match.

`-F` skips invalid prediction records and continues reading the file, retaining
valid records before and after the error. Invalid records do not contribute to
prediction counts. It does not relax experimental-data validation or allow missing prediction files
or missing `F:` file lists.

## Nearest-neighbour context

Use `-N` to append the nearest non-self neighbour distance from a TFDataRecord
file containing serialized `nnbr::NearNeighbours` protos:

```shell
gfp_nearneighbours_single_file_tbb -S dataset.neighbours.tfrecord dataset.gfp
mispredicted -E measured.txt -N dataset.neighbours.tfrecord -d pred1.txt pred2.txt
```

With `-d`, an additional `NN_DIST` column appears in both continuous and
classification output. Normal descriptive output appends `NN_DIST <distance>`.
The continuous `-Y` summary also gains an `NN_DIST` column, and `-T` textprotos
include optional `nearest_neighbour_distance`. Existing ranking, prediction
counts, filtering, and output order are unchanged. Without `-N`, output remains
unchanged.

Identifiers are matched exactly. Missing query IDs and records containing no
non-self neighbours produce `.` in text output; the proto field is absent.
Distinct neighbours with distance zero are retained. Duplicate query IDs are an
error. Neighbour lists are trusted to be sorted nearest first; self-neighbours
are skipped. By default the file contains `NearNeighbours` protos. For indexed
`NearNeighboursIndices` input, specify the format explicitly:

```shell
gfp_nearneighbours_single_file_tbb -S dataset.indexed.tfrecord -S indexed dataset.gfp
mispredicted -E measured.txt -N dataset.indexed.tfrecord -N indexed -d pred1.txt pred2.txt
```

The two `-N` arguments can appear in either order. Specify exactly one filename
and at most one `indexed` qualifier. Format detection is not automatic: choose
the format that was written to the file.

Indexed neighbour IDs are zero-based record positions. Keep the original record
order and all records, including those with no neighbours; do not filter or
reorder an indexed file. Out-of-range retained neighbour indices are errors.
The reader retains only the nearest non-self neighbour per molecule and resolves
its index to a name, including forward references. All output formats behave
the same for equivalent named and indexed input.

With `-v`, report how many records were loaded and coverage among all identifiers
with at least one valid prediction, before prediction-count and error filters:
identifiers with distances, absent IDs, and empty neighbour lists are counted
separately.

The distance describes the reference dataset used to create the supplied file.
For repeated prediction splits, a whole-dataset neighbour file provides structural
context for persistent mispredictions, rather than the exact nearest distance
in each split's training set. A small distance alongside a large prediction error
can suggest a different investigation from a large nearest-neighbour distance;
the tool does not determine the cause.

## Error metrics and selection

Continuous errors default to **mean absolute error**:

```text
mean(abs(measured - predicted))
```

`-w <count>` selects at most that many molecules with the largest mean errors.
`-W <threshold>` retains molecules whose mean error magnitude is at least the
threshold. Use one of these selection options at a time; if both are given,
`-w` takes precedence.

`-s` instead uses **mean signed error**, `mean(measured - predicted)`. Positive
values indicate under-prediction; negative values indicate over-prediction.
With `-s`, both `-w` and `-W` select by the magnitude of that signed mean.
Opposite-sign errors can cancel: errors of +2 and -2 have mean absolute error 2
but mean signed error 0.

`-b` requires `-w` and implies signed errors. It selects from both directions,
approximately balancing accumulated error magnitudes while both signs remain
available, and never exceeds the requested count. One-sided data remains valid.

`-n <count>` requires at least that many **valid prediction records** per
identifier, default 1. Duplicate records within a file count separately: this is
not a count of independent files or models. A filter with no qualifying molecules
succeeds with no data rows. Input containing no valid predictions is an error.

Continuous output sorts by decreasing mean error, then decreasing prediction
count, then ascending identifier. For signed errors this puts positive means
before negative means; selection still considers both directions. `-r` reverses
the final output order, including ties, without changing selected identifiers.

## Output

Default continuous output is a descriptive line per identifier. `-d` produces a
headered table:

```text
ID EXPT AVDIFF
mol1 2 1
mol2 4 0.5
```

`AVDIFF` is the selected mean error metric. `-S molecules.smi` prepends SMILES to
continuous output; the SMILES file contains `smiles id` records without a header.
Missing SMILES are reported and represented by `*`.

`-Y summary.txt` writes a continuous summary with these columns:

| Column | Meaning |
| --- | --- |
| `ID` | Identifier |
| `Expt` | Experimental value |
| `N` | Number of valid prediction records |
| `minpred`, `maxpred`, `avepred` | Prediction range and mean |
| `avediff` | Mean error using the selected metric |
| `maxdiff` | Largest error using the selected metric; in signed mode, the largest signed value |
| `stddiff` | Sample standard deviation of those errors |

A single prediction has `stddiff` 0 by convention. With the default metric this
is spread in absolute errors; with `-s` it is spread in signed errors.

`-T results` writes one single-line `MispredictedContinuous` textproto per
identifier to `results.textproto`. The suffix is appended if needed. Summary and
textproto outputs use the same eligible identifiers and order as primary output.
Output paths that alias input files or each other are rejected before opening.

## Selecting prediction files

Pass filenames explicitly, or use `-R` to match filenames in the current directory:

```shell
mispredicted -E measured.txt -n 5 -d -R '^prediction_[0-9]+\.txt$'
```

Matching is against the filename, not a recursive directory search. Only regular
files (including symlinks to regular files) are considered. The matching filenames
are collected and sorted before output files are created, so a newly created
summary or proto file cannot become an input during the same run. With `-R`,
positional prediction filenames are not used.

Continuous mode also accepts a list file:

```shell
mispredicted -E measured.txt -d F:prediction_files.txt
```

The list contains one filename per line, without a header. Empty lines and lines
starting with `#` are skipped. Paths are relative to the working directory.

## Classification

`-C` compares class-label tokens instead of numeric errors:

```shell
mispredicted -C -E measured_classes.txt -n 2 -d pred1.txt pred2.txt
```

Descriptor output has `ID EXPT N Correct Accuracy`, where accuracy is the fraction
of prediction records matching the experimental class. Default ordering is best
accuracy first, then largest prediction count, class label (length before lexical
order), and identifier. `-r`
reverses that order, putting the least accurate identifiers first.

`-c <cutoff>` enables classification from numeric scores: scores below the cutoff
become `-1`, and scores at or above it become `1`. `-X <label1,label2>` translates
predicted labels to experimental classes by prevalence: the first predicted label
maps to the more common experimental class. This requires exactly two experimental
classes and two distinct predicted labels. Equal class populations are ordered
lexically by their experimental label, so reordering the input cannot reverse
the translation.

Classification does not implement `-w`, `-W`, `-Y`, or `-T`; these options are
rejected before any output file is opened. `-S` supplies the SMILES prefix in
both classification and continuous output, including descriptor rows. An
identifier missing from the supplied SMILES file is written with `*` as its
SMILES and a diagnostic on stderr.

## Support-vector distance options

`-k def` (or `-k default`) reads a distance to the training set immediately
following the selected prediction column. The distance must be finite and in
[0, 1]. It is retained as `NN_DIST` in descriptive, tabular and summary output,
and as `nearest_neighbour_distance` in textproto output. It does not affect
ranking. If an identifier occurs more than once, its last valid record supplies
the distance. This option is only supported for continuous predictions and
cannot be combined with `-N`. Other format names are rejected.

The `-D` distance/error output option is not implemented and is rejected.

Summary and textproto write failures, including final buffer flush failures,
produce a nonzero exit status. Standard deviation is the sample standard
deviation of the selected error metric; constant repeated errors have zero
spread. The shared accumulator uses double-precision running central moments
rather than subtracting nearly equal raw sums of squares.

## Build and regression tests

From `src`:

```shell
bazel build -c opt --config=native //Utilities/General:mispredicted
ruby ../test/run_all_test.rb -bindir "$PWD/bazel-bin/Utilities/General" mispredicted
```

The native configuration is for local builds targeting the host CPU.
