# dicer_fragment_activity

`dicer_fragment_activity` associates algorithmically generated dicer fragments
with activities. It compares molecules that generated a fragment with those
that did not, then optionally writes a JSON report for an interactive static
HTML display. Structures are treated as strings; neither the analyzer nor the
fragment HTML renderer requires LillyMol's Python bindings or RDKit.

## Generate, analyze and display

The activity file needs a header followed by identifier and numeric activity:

```text
ID activity
mol1 5.2
mol2 7.1
```

Generate per-molecule serialized dicer output, analyze it, and build the report:

```shell
dicer -k 3 -I 1 -B serialized_proto -S diced.dat molecules.smi
dicer_fragment_activity -A activity.txt -tfdata -support 10 \
  -json fragment_activity.json diced.dat
dicer_fragment_activity_html.sh fragment_activity.json \
  --top_n 100 -o fragment_activity.html
```

Use dicer's per-molecule `DicedMolecule` output. Aggregate fragment-statistics
files do not contain the molecule memberships required for this analysis.
For text input, generate `dicer -B proto ...` and omit `-tfdata` when analyzing.
Text input has one complete `DicedMolecule` textproto per line. Multiple input
files can be supplied when they contain disjoint molecule identifiers and were
created with consistent dicer settings.

`-json -` writes JSON to stdout, allowing a pipe into the HTML renderer:

```shell
dicer_fragment_activity -A activity.txt -tfdata -support 10 -json - diced.dat |
  dicer_fragment_activity_html.sh - --top_n 100 -o fragment_activity.html
```

The normal build installs `dicer_fragment_activity`. To build it locally, from
`src`, use `bazel build -c opt //Molecule_Tools:dicer_fragment_activity` (add
`--config=native` for local CPU targeting). The HTML scripts are in `contrib/bin`
and require only the Python standard library. The shell wrapper uses `PYTHON`
when set, otherwise `python3`.

## Population and fragment identity

The population contains the unique molecule records actually read, including
molecules with no generated fragments. Activity records not present in the
dicer input are excluded from all baseline statistics and prevalence counts.
Missing activities, repeated molecule identifiers and repeated activity
identifiers are errors. A full molecule name is matched first; if necessary,
its first whitespace-separated token is used to find the activity identifier.
Numeric activities must be finite.

Each fragment contributes at most once per molecule. Repeated fragment entries
are deduplicated, and `DicerFragment.n` is not an activity weight. Fragment
identity is its coordinate-free `usmi` when present, otherwise `smi`. Isotope
labels remain part of the identity. The tool does not parse or recanonicalize
these strings; ordinary dicer output uses unique SMILES by default.

An absent fragment means it was not generated under the chosen cut, size,
stereochemistry and isotope settings. This does not necessarily mean that its
unlabelled substructure is absent from the parent. Fragment-subset redundancy
is deliberately left untouched.

Corrupt serialized protos, checksum errors and truncated records cause a
nonzero exit status. Input failures do not produce a partial report.

The shared TFDataRecord reader now distinguishes truncated frames from clean
EOF and reads packed scalar fields without alignment assumptions. The writer
also includes the required data CRC for empty payloads. Old zero-length records
written without that CRC were nonstandard and must be regenerated; normal
per-molecule dicer records have nonempty payloads.

## Analyzer options

| Option | Meaning |
|---|---|
| `-A file` | Headered whitespace-separated activity file. |
| `-tfdata` | Read serialized TFDataRecord protos instead of text. |
| `-support n` | Require at least `n` matched molecules; default retains all observed fragments. |
| `-C n` | Maximum fragment atom count, using the proto's `nat` field. |
| `-json file` | Write JSON instead of the traditional table; `-` means stdout. |
| `-unmatched_support n` | Require at least `n` unmatched molecules in JSON; default 1. |
| `-bins n` | Number of common histogram bins; default 20. Constant activities use one bin. |
| `-examples n` | Maximum parent examples per JSON feature; default 3, or 0 for none. |
| `-v` | Report molecule, fragment and unused-activity counts. |

JSON contains all qualifying fragments, matched/unmatched counts, minimum,
maximum, mean, median and quartiles, histograms on a common activity scale,
and signed matched-minus-unmatched mean and median differences. Quantiles use
linear interpolation at `(n - 1) * percentile`. Parent examples contain SMILES,
identifiers and activities; they are not atom-mapped or highlighted.

Exact statistics use compact molecule-membership lists. Matched/absent groups
are accessed through global activity ranks: absent quartiles do not require
scanning every molecule for every fragment, and absent histograms are obtained
by subtraction from the global histogram. Each feature is written to JSON
individually. JSON mode uses additional membership memory; table mode does not
retain memberships or parent SMILES.

Without `-json`, the existing columns remain:

```text
Frag Parent Natoms N Fraction Minval Ave Maxval Diff
```

`Fraction` uses the processed population. The table's `Diff` retains its meaning
of matched mean minus population mean; JSON's `delta_mean` compares the matched
and unmatched groups. Both output forms have deterministic fragment-SMILES
ordering. Output is rejected if its path would overwrite an input file.

## Ranking and HTML limits

The HTML renderer sorts by absolute median difference by default, retaining the
100 most discriminating fragments. Use `--top_n 0` to show all eligible rows.
JSON is not reduced by this display limit, so alternative views can be generated
without rereading dicer data.

| Option | Meaning |
|---|---|
| `--top_n n` | Maximum HTML rows after filtering/ranking; default 100, 0 means all. |
| `--sort abs_median` | Rank separation in either direction; the default. |
| `--sort abs_mean` | Rank absolute mean differences. |
| `--sort median` | Rank signed median differences, largest first. |
| `--sort mean` | Rank signed mean differences, largest first. |
| `--sort support` | Rank by matched molecule count. |
| `--ascending` | Reverse score direction, e.g. most negative signed differences first. |
| `--min_matches n` | Additional matched-group minimum; default 1. |
| `--min_unmatched n` | Additional unmatched-group minimum; default 1. |
| `--title text` | Report title. |
| `-o file` | Required HTML output; `-` means stdout. |

Ties are resolved by fragment SMILES and feature name. Additional HTML filters
cannot restore features omitted by the analyzer's support or size limits. The
renderer reads the complete JSON; `--top_n` limits rendering, not JSON parsing
memory. The report uses the same HTML rendering code as
`substructure_exploration.py`, including text filtering, sortable columns and
matched/unmatched activity histograms.

These are exploratory associations within the analyzed dataset; correlated
fragments may have identical distributions.
