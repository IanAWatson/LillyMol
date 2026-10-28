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

The default catalogue is `contrib/data/substructure_exploration.textproto`. Its
83 features include functional groups, halogens, common heterocycles, ring-system
topology, aromatic linkers and potentially reactive groups. Definitions overlap:
a lactam is also an amide, and a carbamate may also match ester and amide queries.
These are descriptive features, rather than disjoint chemical classifications.

Some definitions deserve particular attention:

- `cf3` requires exactly three fluorines attached to a carbon, with a fourth
  non-fluorine neighbour. It excludes CF4.
- `heavy_halogen` matches Cl, Br or I.
- `methoxy` matches `[CD1H3]-[OD2]-[c,CX4]`: methyl ethers attached
  to aromatic or saturated aliphatic carbon, excluding methyl esters and
  non-aromatic unsaturated carbon attachments.
- `fluorine_not_cf3` matches individual fluorines outside CF3, including aryl
  fluorine, monofluoroalkyl and difluoroalkyl fluorines.
- The primary, secondary, tertiary and cyclic amine queries use separate,
  commented `environment_no_match` messages to exclude common amide-like,
  amidine/guanidine and directly heteroatom-linked nitrogens. Aniline is also
  available as a separate, overlapping feature.
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
chemically distinct groups for every complex motif. The six inline ring-system
predicates are **presence only**: they report zero or one match even when several
systems qualify, and highlight one atom belonging to a matching system. This
avoids counting all atoms of a qualifying system as separate occurrences.

## Custom queries and reports

Use `--config queries.textproto` to replace the catalogue. Each entry needs a
unique name and one of SMARTS, a query filename, or an inline full LillyMol
`SubstructureQuery` proto:

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
are resolved from the working directory. Legacy `ncon` and `charge_assigner`
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

The molecule-independent report renderer is shared with
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
