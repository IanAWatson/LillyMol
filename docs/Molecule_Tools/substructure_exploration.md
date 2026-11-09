# substructure_exploration

`substructure_exploration.sh` compares activity values for molecules that contain
or lack medicinal chemistry features. It writes JSON with summary statistics,
activity distributions, match counts and labelled examples. It can also produce
a static HTML report.

```shell
substructure_exploration.sh analyse -A activity.txt molecules.smi \
  -o exploration.json --html exploration.html
```

The activity file must have a header, followed by identifier and activity columns:

```text
ID activity
mol1 6.2
mol2 7.1
```

Molecule identifiers are the first token of each SMILES name. Activity files may
be whitespace separated, or comma separated when the filename ends in `.csv`.
Use `--okmissing_activity` to skip molecules without activity values.

The wrapper locates the Python script and default catalogue independently of the
current directory. It uses `$LILLYMOL_HOME/run_python.sh`; if `LILLYMOL_HOME` is
unset, it infers the repository root from the wrapper location. Use a Python
interpreter matching the staged LillyMol bindings, with the `protobuf` package
installed. Set `PYTHON=/path/to/python` when needed.

## Matching and standardisation

All structure processing and substructure searching use LillyMol. By default,
all LillyMol chemical standardisations are enabled before matching. This makes
common neutral and ionised representations of acids and amines comparable;
it does not assign physiological charge states. Permanent charges, such as
quaternary ammonium, remain charged. Examples in the report show the standardised
structures. Use `--no-standardize` to search the input structures directly.

The default selection file is `contrib/data/substructure_exploration/queries`.
It lists 85 standalone LillyMol query files in the same directory. The features include functional groups, halogens, common heterocycles, ring-system
topology, aromatic linkers and potentially reactive groups. Definitions overlap:
a lactam is also an amide, and a carbamate may also match the ester query.
These are descriptive features, rather than disjoint chemical classifications.

Some definitions deserve particular attention:

- `cf3` matches `F-[CD4](-F)-F`. It also matches CF4; unique matched atom
  sets give four matches for CF4, one for each choice of three fluorines.
- `heavy_halogen` matches Cl, Br or I.
- `methoxy` matches `[CD1H3]-[OD2]-[c,CG0]`: methyl ethers attached
  to aromatic or saturated aliphatic carbon, excluding methyl esters and
  non-aromatic unsaturated carbon attachments.
- `fluorine_not_cf3` matches fluorines bonded to aliphatic carbon outside CF3,
  including monofluoroalkyl and difluoroalkyl fluorines. Aromatic fluorine is
  covered separately by `aryl_fluoride`. The non-CF3 query also matches the
  fluorines in CF4, so it overlaps the `cf3` feature for CF4.
- `amide` uses `[CX3D3T2](=O)-[NX3]`, excluding formamides, ureas and
  carbamates from this general amide feature. The separately named amide
  subtypes can overlap other carbonyl features.
- `ether` and `ether_cyclic` require two saturated carbons around the oxygen;
  they distinguish non-ring and ring oxygen respectively. Aryl ethers have
  their own feature, and ester/lactone oxygens do not satisfy these queries.
- The primary, secondary, tertiary and cyclic amine queries use separate,
  commented `environment_no_match` messages to exclude common amide-like,
  amidine/guanidine and directly heteroatom-linked nitrogens. Aniline is also
  available as a separate, overlapping feature.
- `sulphur_carbon_carbon` matches `[#6]-[SX2]-[#6]`, allowing aromatic as
  well as aliphatic carbon neighbours.
- `saturated_furan` identifies an unfused saturated five-membered oxygen ring.
- `thiazol_2_one` identifies the carbonyl motif on an unfused five-membered
  aromatic ring with adjacent sulfur and nitrogen atoms.
- Common named ring queries use `/IWfss1` to restrict them to unfused rings.
- Terminal, internal and standalone aromatic systems have respectively one,
  at least two, or zero non-spinach connections. A pendant acyclic substituent
  does not make an otherwise standalone aromatic system terminal.
- `fused_aromatic_aromatic_2` requires exactly two aromatic rings in the system;
  `fused_aromatic_aliphatic_2` requires one aromatic and one non-aromatic ring.
- `biaryl_bond` explicitly excludes ring bonds, including Kekule single bonds
  within aromatic systems.

Matching counts use unique matched atom sets to suppress symmetry-equivalent
embeddings. They are counts of query matches, which need not equal counts of
chemically distinct groups for every complex motif. The six ring-system
predicates are **presence only**: they report zero or one match even when several
systems qualify, and highlight one atom belonging to a matching system. This
avoids counting all atoms of a qualifying system as separate occurrences.

## Custom queries and reports

### Selection files and standalone queries

Use `--config <selection>` to choose a subset. A selection file lists one query
pathname per line, in the desired report order. Blank lines and lines beginning
with `#` are ignored. Each listed file is a normal LillyMol `SubstructureQuery`
textproto with its own name, for example `acid.textproto`:

```textproto
name: "acid"
query {
  smarts: "[OD1H]-C=O"
}
```

Paths are relative to the selection file, not the current working directory.
Absolute paths, spaces in filenames, `~`, and environment variables such as
`${LILLYMOL_HOME}` are supported. Each line is a complete pathname, without shell
quoting or inline comments. Query names must be unique across the selection.

For example, a file called `my_queries` could contain:

```text
# A small feature subset
${LILLYMOL_HOME}/contrib/data/substructure_exploration/acid.textproto
${LILLYMOL_HOME}/contrib/data/substructure_exploration/primary_amine.textproto
```

Run it with:

```shell
substructure_exploration.sh analyse --config my_queries -A activity.txt molecules.smi
```

The supplied `functional_groups.txt` and `rings.txt` selections provide starting
subsets. Copy a selection and edit its entries to experiment, or edit an
individual query file without changing every selection that uses it. Standalone
files also work with other LillyMol tools:

```shell
tsubstructure -q "PROTO:$LILLYMOL_HOME/contrib/data/substructure_exploration/acid.textproto" molecules.smi
```

### Existing textproto configurations

The `SubstructureExploration` textproto format remains supported via
`--config queries.textproto`. The old
`contrib/data/substructure_exploration.textproto` path is retained as a small
compatibility configuration referencing the standalone query files; definitions
are maintained only in those individual files.

Each textproto entry can specify SMARTS, a query filename, or an inline full
LillyMol `SubstructureQuery` proto. An explicit name overrides the query file's
name. When a filename entry omits the name, it is read from that file's textproto.
Relative filenames and environment variables are resolved as in selection files.
For legacy MSI query files, continue to supply an explicit name.


```textproto
query {
  name: "CF3"
  smarts: "[CD4;!$(C(F)(F)(F)F)](F)(F)F"
}
query {
  name: "terminal_aromatic_system"
  query {
    query {
      ring_system_specifier {
        base { set_global_id: 1 }
        number_non_spinach_groups: 1
        min_aromatic_ring_count: 1
        non_aromatic_ring_count: 0
      }
      smarts: "[/IWgid1]"
      max_matches_to_find: 1
    }
  }
}
```

For a query file, use `fname: "path/to/query.textproto"`. Relative query filenames
are resolved from the configuration file directory. Legacy `ncon` and `charge_assigner`
metadata are recorded but not applied; the script warns if they are supplied.

`--min_matches 1` is the default: only queries matching at least one molecule
with an activity value are included in the JSON and HTML output. Increase the
threshold to require more support, or use `--min_matches 0` to retain every
query, including those with no matches. This counts molecules, not embeddings
or retained display examples. The global activity statistics still use all
molecules with activity. Metadata records the total queries, queries reported
and threshold. HTML rebuilt from JSON uses the queries retained in that JSON.

`--max_examples 3` limits examples retained per feature. `--bins 20` controls the
common activity histogram bins, and `--count-buckets '0,1,2,3+'` controls the
multiplicity profiles. Rebuild an HTML report from existing JSON with:

```shell
substructure_exploration.sh html exploration.json -o exploration.html
```

When rendering existing JSON, `--sort abs_median --top_n 100` selects the
largest activity separations in either direction. Other sort keys are
`abs_mean`, `median`, `mean` and `support`; `--ascending` reverses score direction.
`--min_matches` and `--min_unmatched` can further filter the retained queries.
By default, the `html` command preserves the stored order and includes all
stored queries. These options change HTML only and cannot restore queries
filtered out during analysis.

### Pagination

HTML reports display 50 query rows per page by default, with Previous/Next
buttons and a page count. `--page_size n` sets the number of visible rows;
`--page_size 0` shows every retained row at once. This option is available for
both `analyse --html` and the `html` command:

```shell
substructure_exploration.sh html exploration.json --page_size 25 -o report.html
```

Sorting and filtering operate on all retained rows before pagination and reset
the view to page one. Pagination changes browsing only: all retained rows and
images remain embedded in a single offline HTML file. Use `--top_n` when
rendering existing JSON to limit total output size and depiction-generation
work, independently of `--page_size`.

### Hide median statistics

Use `--nomedian` to suppress median activity values and the Delta Median column
in HTML reports. It also hides the global median and median sorting controls;
the global summary, sorting control and multiplicity profiles use means instead.
The JSON still contains all computed median statistics, so it can be rendered
again with medians later.

```shell
substructure_exploration.sh analyse --nomedian -A activity.txt --html report.html molecules.smi
substructure_exploration.sh html exploration.json --nomedian -o report.html
```

With `--nomedian`, ranking requested through `--top_n` defaults to `abs_mean`
instead of `abs_median`. An explicit `--sort` choice is still respected. Without
a ranking request, the stored row order is preserved as usual.

### Structure examples and long queries

Matched examples are shown as embedded RDKit SVG depictions, with matched atoms
highlighted in green. Isotope-1 markers from the stored labelled examples are
hidden in the drawing; the stored SMILES are unchanged and can be expanded
under each image for copying. Identifiers, activities and match counts remain
visible beneath the structures.

To enable depictions, install RDKit in the Python environment used by the tool:

```shell
uv pip install rdkit
# Alternatively:
python -m pip install rdkit
```

RDKit is loaded only for HTML structure rendering. LillyMol still performs
standardisation and substructure searching, and analysis without HTML does not
require RDKit. If RDKit is unavailable, or cannot depict a particular LillyMol
SMILES, the example displays a message and its SMILES instead. Generated HTML
is self-contained, including its images, and can be viewed offline.

Query/source definitions longer than 100 characters show a shortened preview
and an expandable **Full query** section. Text filtering still searches the
complete definition, including the hidden part. JSON retains the full source.

The report renderer is shared with
[dicer_fragment_activity](dicer_fragment_activity.md), which analyzes generated
fragments rather than dictionary queries.

The report compares groups within this dataset. Differences may reflect
correlated features and chemical series as well as the queried feature.

## Development

The generated `contrib/bin/substructure_exploration_pb2.py` accompanies the
script. After changing its schema, regenerate it from `src` with a `protoc`
compatible with your Python protobuf runtime:

```shell
protoc -I ../contrib/bin -I . --python_out=../contrib/bin \
  ../contrib/bin/substructure_exploration.proto
```

The imported `Molecule_Lib` Python protos are supplied by the LillyMol Python
build/staging process. Run the catalogue and wrapper regression tests with:

```shell
$LILLYMOL_HOME/run_python.sh \
  $LILLYMOL_HOME/contrib/bin/test/substructure_exploration_test.py
```

Each catalogue entry has a matching and nonmatching example in
`contrib/bin/test/substructure_exploration/examples.json`.
