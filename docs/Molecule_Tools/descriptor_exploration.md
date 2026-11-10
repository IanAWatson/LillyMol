# descriptor_exploration

`descriptor_exploration` compares measured outcomes with an `iwdescr` table.
It ranks properties for inspection, shows linear and nonlinear relationships,
and distinguishes an interesting-looking winner from one supported by held-out
predictions and shuffled-activity comparisons. It writes a complete JSON
analysis and an interactive, self-contained HTML report.

## Quick start

```shell
uv pip install numpy scipy scikit-learn
iwdescr.sh molecules.smi > molecules.w
descriptor_exploration.sh analyse molecules.w -A activity.txt \
  -o exploration.json --html exploration.html
```

The app itself uses Python with NumPy, SciPy and scikit-learn (1.2 or newer).
It does not require LillyMol Python bindings when analysing existing tables.
The shell wrapper uses the checked `run_python.sh` environment when staged
bindings are installed, and otherwise uses `${PYTHON:-python3}`. Set `PYTHON` to
an interpreter in your environment as with other LillyMol scripts.

Both inputs require headers. The first column is always a string identifier;
records are joined by identifier, never by position. The measurement defaults
to the second activity column. Use `--activity_column assay_name` to select
another column by header name.

```text
Name w_natoms w_xlogp w_nvrtspsa
mol1 20 1.2 45
mol2 27 3.5 62
```

```text
ID activity
mol1 5.2
mol2 7.1
```

These two records illustrate the format; meaningful analysis requires a larger
dataset. Whitespace-delimited files, `.csv`, `.tsv`, and their gzip-compressed
forms are supported. Descriptor columns must be numeric. Missing descriptor
values can be `.`, `NA`, `NaN`, `nan`, `null`, `None`, or an empty CSV/TSV field.
Duplicate identifiers, ragged tables, non-finite numerical values and censored
measurements such as `>10` are rejected. Censored values must be handled
explicitly before analysis; replacing them with the limit assumes an exact
measurement and changes the interpretation.

Missing measurements cause an error unless `--allow_missing_activity` is given.
Extra activity identifiers are counted but ignored. The JSON records matched,
skipped and unused counts. All statistical calculations use the matched set.

## Activity scale and intended direction

```shell
descriptor_exploration.sh analyse molecules.w -A measurements.csv \
  --activity_column IC50_molar --transform neglog10 --direction higher \
  --units pIC50 -o exploration.json --html exploration.html
```

Transforms are `none` (default), `log10`, and `neglog10`. Log transforms require
positive measurements. No concentration-unit conversion is implicit: pIC50
requires molar input. `--direction higher|lower` describes the desirable
**transformed** response and is displayed in the report. It does not change
association strengths or turn an observed relationship into a causal effect.
`--units` labels the transformed measurement scale.

## Descriptor filtering and buckets

Constant descriptors are omitted. Default near-constant filters omit features
whose most frequent value accounts for at least 99% of observed records, or
whose standard deviation is at most `1e-12 * max(abs(value))`. These are
configurable using `--max_mode_fraction` and `--min_relative_sd`. A rare value
can still be scientifically interesting: adjust the filters deliberately for
such analyses. The relative-spread threshold is a numerical guard, not a
universal measure of chemical relevance. Omitted descriptors and reasons
remain in JSON and HTML. LillyMol's `notenoughvariance` can also be used to
prepare a descriptor table before running this app.

At least eight observed measurements per descriptor are required by default
(`--min_samples`, minimum permitted 4). Features with more than 20% missing
values are omitted (`--max_missing_fraction`). Univariate analyses exclude
missing values feature by feature; no activity-based imputation is performed.

Requested bucket counts default to `--buckets 2,4,8`. Identical descriptor
values are never split. Features with few distinct values initially use those
values as groups; other features use approximately equal-population groups.
Adjacent groups smaller than `--min_bucket_size 5` are merged, preferring the
less populous neighbor. If this leaves just one group, its descriptive η² and
bucket prediction gain are zero. Requested and actual bucket counts are shown.
These rules apply to all features, including integer counts.

Each group shows population, observed descriptor range, mean activity, optional
median, activity standard deviation, and an approximate 95% t interval for the
mean. Those intervals assume independent observations and do not account for
chemical-series dependence or simultaneous comparison of many buckets.
`--nomedian` hides medians in HTML but retains them in JSON.

## Association measures and nonlinear fits

For each retained descriptor the app computes:

- Signed Pearson and Spearman correlations.
- Mutual information (MI) against activity. Low-cardinality integer descriptors
  use the discrete-feature estimator when every level has at least two records;
  other features use the continuous estimator. The choice is recorded in JSON.
- Bucket η²: variance explained by differences between group means. It is
  descriptive and grows optimistically as groups become more flexible.
- Held-out prediction gains for a straight-line fit, a regularized spline fit,
  and each requested bucket resolution.
- A descriptive full-data spline curve, its activity range, and the fitted
  activity difference between descriptor quartiles. These effect magnitudes
  use transformed activity units, rather than a unitless correlation scale.

Splines join local polynomial pieces smoothly. They can represent a middle
optimum or a threshold even when signed correlations are close to zero.
Up to five quantile knots are used, duplicate knots are merged, and spline
degree is reduced for features with fewer distinct values. Ridge
regularization uses a fixed penalty of 1. No knot-count or penalty search is
performed, so the held-out observations are not used to tune these settings.

Repeated shuffled folds default to `--folds 5 --repeats 3`. Folds are reduced
for small observed subsets so each test fold has at least two records.
For each fold, spline knots, scaling, bucket boundaries, and bucket means are
learned on training molecules only. The prediction gain is:

```text
1 - sum(held-out prediction squared errors)
    / sum(held-out squared errors from the training-fold mean)
```

Zero means no improvement over the mean; negative is worse. Gains are computed
from pooled errors within each repeat and then averaged. The repeat values
are shown to make instability visible. The displayed red spline is fitted to
all observations for description: it is **not** the held-out prediction curve.
Gray points may be reduced by `--max_points 400` for display only. Bucket
statistics, scores, correlations and molecule exports always use all relevant
observations. Hover over a plotted point for its identifier and values.

## Composite score and evidence

The composite inspection score is a bounded heuristic, not a probability:

```text
trend       = max(Pearson², Spearman²)
dependence  = ((1 - exp(-2*MI)) + max_bucket_eta²) / 2
score       = 100 * (0.2*trend + 0.2*dependence + 0.6*clip(spline_gain, 0, 1))
```

The MI transformation supplies a bounded scale; it is not a universal conversion
of MI into explained variance. Related measures share a component so that
several versions of the same trend do not each get a full vote. Negative
prediction gains remain visible in the table even though their composite
prediction contribution is zero. Forest importance does not enter this score.

`--permutations 19` shuffles activities while retaining the descriptor matrix.
The same shuffle is used across descriptors. Each shuffled response is scored
with the same correlations, MI, buckets and validation procedure. A descriptor's
permutation p compares against its own shuffled scores. The displayed
**familywise p** compares against the maximum score across all retained
descriptors in each shuffle. Both use `(1 + exceedances)/(1 + permutations)`.
This helps distinguish a real-looking winner from the best accident among many
features. Fixed `--seed 2026` makes shuffles and splits reproducible.

With 19 shuffles the finest p resolution is 0.05. Use `--permutations 99` or more
for a finer comparison, at greater computational cost. `--permutations 0`
disables calibration and prevents a supported label. Shuffle calibration
assumes exchangeable measurements and is exploratory, not independent
confirmation. Correlated series, batches, and repeated measurements can violate
that assumption. Row-wise folds can also overstate transfer to new chemistry.
Series-aware analysis is intentionally deferred in this first version.

Evidence labels are deliberately separate from ranking:

- **Supported:** familywise p ≤0.05; mean spline gain at least
  `--min_cv_gain 0.02`; positive gain in every repeat; descriptive curve range
  at least `--min_effect`; and at least 30 observed records.
- **Tentative:** gain, repeat stability and effect checks pass, but at least
  one of the stronger evidence checks does not.
- **Weak:** insufficient held-out improvement, stability or effect size.

These thresholds are inspection rules, not guarantees. `--min_effect` defaults
to zero because practical relevance depends on the measurement. Set it in
transformed activity units to require a scientifically useful difference.
Small-sample warnings remain visible even if a descriptor ranks first. A weak
winner is explicitly identified in the report summary.

## Correlated properties and the random forest

Descriptor correlations are learned from the current matched dataset.
Pairwise-complete Spearman correlations with absolute value at least
`--correlation_threshold 0.85` are listed with each descriptor. These are
neighbors above a threshold, not claims of independent property effects or
transitive clusters. Their sample counts are recorded.

A multivariable random forest runs by default (`--no_forest` disables it).
It uses `--trees 100`, a minimum leaf size of 3, the same repeat/fold settings,
and train-fold median imputation. Feature importance is the increase in
held-out mean squared error after permutation, normalized by the fold's
training-mean baseline error, averaged over folds/repeats. It can be negative.
The forest's overall held-out prediction gain is displayed; a model that does
not improve on the baseline gets an explicit warning.

Correlated descriptors can substitute for one another, making each look
unimportant. Use the correlation neighbors alongside this result. The forest
column asks about contribution in a multivariable model; univariate scores ask
about association in isolation. Neither establishes what will happen if a
chemist changes the corresponding property.

`--jobs` controls forest and forest-permutation parallelism. A single worker is
the default. The optional sd of forest importance in JSON describes variability
across fold/repeat means, not a confidence interval.

## Molecules, exports and report size

```shell
uv pip install rdkit   # optional depictions only
# Run with staged LillyMol bindings and a matching Python ABI:
$LILLYMOL_HOME/run_python.sh $LILLYMOL_HOME/contrib/bin/descriptor_exploration.py \
  analyse molecules.w -A activity.txt --smiles molecules.smi \
  -o exploration.json --html exploration.html
```

Structure enrichment uses LillyMol to read SMILES and determine atom counts.
RDKit is used only for optional depiction. The SMILES file has SMILES followed
by identifier (first name token); missing examples remain identified by name.
No structure is standardized, transformed or used to alter the descriptor
values by this app. Generate consistent descriptors and SMILES upstream.

Select a bucket's **Inspect / export** button to see examples near the lowest,
midpoint and highest activities. Ties prefer smaller structures when supplied.
Export includes every member identifier and its transformed activity, including
molecules omitted from the point display. Images are limited globally by
`--max_structures 100`; other available examples show SMILES. The HTML works
offline without external JavaScript, image services or plotting libraries.

JSON retains all analysed descriptors. HTML defaults to the highest-scoring
50 (`--top_n 0` includes all) with 20 descriptors per page (`--page_size`).
Filtering and sorting act on all embedded descriptors before pagination.
The **Export displayed ranking** button exports all rows matching the current
filter, not just the current page. Pagination does not reduce the bytes already
embedded; use `--top_n`, `--max_points`, and `--max_structures` to bound rendering.
Bucket memberships can still be large, since complete molecule exports are
retained.

Rebuild HTML without repeating analysis:

```shell
descriptor_exploration.sh html exploration.json -o alternate.html \
  --top_n 100 --page_size 25 --nomedian
```

Input/output path collisions are rejected. Previously generated JSON records
all filtering and calibration settings used for that analysis. The app never
silently substitutes RDKit for LillyMol chemistry or treats censored values as
exact observations.

## Development checks

```shell
python3 contrib/bin/test/descriptor_exploration_test.py
```

Run from the repository root using the statistical Python environment. Tests
cover linear and U-shaped signals, noise-only winners, sparse integer groups,
missing values, deterministic shuffle calibration, train-only preprocessing,
forest importance, input joins, output safety, and HTML generation.
