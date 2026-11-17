# descriptor_ranking

`descriptor_ranking.sh` shows which computed molecular descriptors vary with
activity in an SAR dataset. It reads an `iwdescr` table and an activity file and
writes one self-contained HTML report that can be emailed as it is. The report
ranks the descriptors and plots how activity changes across 8 buckets of each
one. It is meant for a quick, qualitative look, and runs in seconds to a minute
or so for a few thousand molecules.

It is one of three tools for exploring an SAR series; see
[SAR series exploration](../Workflows/sar_series_exploration.md) for how they fit together.

```shell
iwdescr.sh molecules.smi > molecules.w
descriptor_ranking.sh -A activity.txt -o descriptors.html molecules.w
```

The activity file has a header, then an identifier and an activity column:

```text
ID activity
mol1 6.2
mol2 7.1
```

Both files are joined by identifier (never by position); molecules that are in
only one file are ignored. The identifier is the first token of the SMILES name,
which is also what `iwdescr` writes in its `Name` column. Files are whitespace
separated, or comma / tab separated if the name ends `.csv` / `.tsv`. The
measurement is the second column of the activity file; use `--activity_column`
to pick another by header name. For potency data such as IC50, use
`--transform neglog10` (or `log10`) so the buckets are on a sensible scale.

## Requirements

The script uses Python with numpy, pandas, scipy, scikit-learn and matplotlib. It does
not use the LillyMol python bindings. The wrapper runs `python3`; set
`PYTHON=/path/to/python` to use an interpreter that has these packages:

```shell
uv pip install numpy pandas scipy scikit-learn matplotlib
```

## What the report shows

For each descriptor:

1. The values are cut into quantile buckets (`--buckets`, default 8). Tied values
   are never split, so a descriptor with few distinct values has fewer buckets.
2. **η²** is the fraction of the activity variance explained by the bucket means,
   from 0 to 1. Unlike a correlation it also detects non-monotonic effects, such as
   activity peaking at intermediate logP.
3. A single random forest is fitted to all descriptors. **RF** is the descriptor's
   impurity importance, its share of the forest's split credit. The importances for all
   descriptors sum to 1. The report also gives the forest's out-of-bag R²; if that is
   near zero, a warning says the importances are weak evidence.
4. **ρ** is the Spearman correlation with activity. Its sign gives the direction of the
   trend, but a U-shaped effect can give a value near zero.
5. The **score** is the mean of the η² and RF percentile ranks (1.0 means best on both).
   It orders the report; it is not a statistical measure.

Beside η² and RF the ranking table shows their ranks, e.g. `0.485 (#1)`. The raw value
tells you how far apart two descriptors are, and the rank tells you where a value sits.
When the two ranks differ a lot, the bucket view and the forest disagree. Hover over a
column heading for its definition.

Each plot shows the bucket means with 95% confidence intervals (red), over a sample of
the molecules (grey), and the dashed line is the overall mean. The axis labels give each
bucket's range of descriptor values and its molecule count. The descriptor's description
from [iwdescr](iwdescr.md) is shown above each plot.

## Correlated descriptors

`iwdescr` descriptors are highly correlated: if molecular weight drives activity, the atom
count will too, and a second plot would add nothing. Descriptors with a Spearman
correlation of at least `--cluster_rho` (default 0.9) are clustered, using average
linkage.

- Only the best-scoring member of each cluster is plotted. The others are listed beneath
  the plot with their correlation to it: `also: w_amw (+0.98), w_nrings (+0.91)`.
- RF importance is summed over the cluster, because correlated descriptors share
  credit and a dominant property would otherwise look modest. η² is the best value in the cluster.
- `--prefer natoms,amw` makes the first of those the representative when it scores within
  0.05 (in percentile rank) of the best member, so a more interpretable descriptor can be shown.
- Every descriptor stays in the "All descriptors" table with its cluster number and role.

A highly correlated descriptor can still carry information the representative lacks, so a
collapsed descriptor is tested for it. The effect of the representative is removed from
both the activity and the descriptor, and η² is computed on what is left (the *residual
η²*). If it exceeds `--residual_eta` (default 0.05), the descriptor is **promoted**: it gets
its own plot, marked "adds to *representative*". The threshold is raised automatically for
small datasets, to `3 × (buckets − 1) / (n − 1)` if that is larger.

Promotion is a screen, not a test. On a dataset where everything correlates with size,
expect several promoted descriptors; raise `--residual_eta` if there are too many plots.

## Output

`-o report.html` is required. It has the activity histogram, the ranking table, one plot per
cluster (the top `--top`, default 25) and, collapsed at the end, the table of all
descriptors and the descriptors dropped. `--csv all.csv` also writes every descriptor:

| Column | Meaning |
|---|---|
| `descriptor` | Descriptor name. |
| `eta2`, `rf_importance`, `spearman_with_activity` | As above, for this descriptor alone. |
| `score` | Percentile-rank score among all descriptors. |
| `cluster` | Correlation cluster number. |
| `representative`, `promoted` | Whether it was plotted for its cluster, or promoted. |
| `residual_eta2` | For collapsed descriptors, the residual η² described above. |

The `score` in the CSV ranks individual descriptors. The `score` in the report ranks the
plots (clusters), with a cluster's summed RF, so the two differ.

## Options

| Option | Meaning |
|---|---|
| `-A file` | Activity file (required). |
| `-o file` | HTML report (required). |
| `--csv file` | Also write all descriptors and their statistics. |
| `--activity_column name` | Activity column, default the second column. |
| `--transform t` | `none` (default), `log10` or `neglog10`. The log forms need positive activities. |
| `--buckets n` | Quantile buckets per descriptor, default 8. |
| `--trees n` | Random forest trees, default 300. |
| `--top n` | Number of plots to draw, default 25. |
| `--cluster_rho r` | Spearman \|ρ\| at which descriptors cluster, default 0.9. |
| `--residual_eta e` | Residual η² needed for promotion, default 0.05. |
| `--prefer a,b` | Descriptors to favour as representatives when scores are close. |
| `--max_missing f` | Drop descriptors missing in more than this fraction of molecules, default 0.2. |
| `--max_mode f` | Drop descriptors where one value has this fraction of molecules, default 0.99. |
| `--jobs n` | Threads for the forest, default all cores. |
| `--seed n` | Random seed, default 1; results are reproducible for a given seed. |
| `--title text` | Report title. |

Descriptors dropped as constant or mostly missing are listed in the report. Missing values
(`.` in `iwdescr` output) in the remaining descriptors are replaced by the median.

## Caveats

- This is exploratory. The ranking tells you what to look at, not what causes activity.
- Descriptors in an SAR series are confounded with the chemical series itself, and
  with each other. A descriptor can rank highly because it marks a sub-series.
- Impurity importance favours descriptors with many distinct values, and it is calculated on
  the training data, so check the out-of-bag R².
- The bucket plots use all molecules, and there is no held-out validation. The
  slower and more elaborate [descriptor_exploration](descriptor_exploration.md) adds
  cross-validation and permutation calibration if you need that.
