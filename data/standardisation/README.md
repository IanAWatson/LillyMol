# Optional external chemical standardisations

These rules are opt-in additions to LillyMol chemical standardisation. They are
not enabled by `-g all` alone. To apply the built-in rules and all distributed
specifications together, use the convenience wrapper:

```shell
standardise_all.sh -o smi -S standardised input.smi
```

[`standardise_all.sh`](../../contrib/bin/standardise_all.sh) runs `fileconv`
with `-g all`, then the guanine and hypoxanthine `EXT:` specifications, and
passes through all user arguments. For example, add `-l` for largest-fragment
processing or use `-o usmi -S -` to write unique SMILES to stdout. Its external
targets apply after built-in processing; the scope restrictions below remain.

To enable an individual specification directly:

```shell
fileconv -g all -g EXT:${LILLYMOL_HOME}/data/standardisation/guanine.textproto \
  -o smi -S standardised input.smi
```

Omit `-g all` to apply only the external transformations. Query filenames in the
specification are resolved relative to that file, so keep the three guanine
files together. Each Standardisation proto occupies one line. Comments and atom
mapping explanations are in the query files. The target SMILES atom order must
agree with the query atom order; it is not an independently mapped reaction.

## Guanine

The templates choose the neutral amino/keto form of a guanine core. They handle
supported hydroxy/oxo and amino/imino forms, and 7H/9H proton placement. The
entire fused purine core is specified, including both exocyclic heteroatoms.
The default target has the hydrogen on N9 when that nitrogen is unsubstituted.
For N7-substituted compounds, a separate template retains that substitution
and imposes the compatible bond pattern. N9 and N1 substituents are retained.
Isotope labels and stereochemistry outside the transformed core are retained.

These first rules deliberately require:

- Neutral atoms throughout the matched core, enforced at runtime by
  `require_neutral_matched_atoms: true` in each query.
- Exactly two fused rings, with an aromatic purine core.
- A terminal exocyclic nitrogen (no N2-substituted amino groups).
- An unsubstituted C8 (two heavy-atom connections).
- A terminal oxygen, not an ether or sulfur analogue.

They do not move explicit hydrogen atoms. Use implicit hydrogens at sites
involved in tautomerisation. More extensively fused, charged, partially
hydrogenated, C8-substituted and N2-substituted analogues are outside this
initial scope. Other purines such as adenine, hypoxanthine, xanthine and
isoguanine are not targets. These restrictions can be extended after separate
validation; they should not be relaxed merely to increase match counts.

Restrictions apply to the molecule as it reaches the external rules;
`-g all` can neutralise a charged input first and make it eligible.

External transformations run after built-in rules. In particular, built-in
imidazole processing can choose the other N-H placement; the external guanine
template then imposes the selected final form. This intentionally broadens
structural equivalence. External templates override the built-in `ENOL`
preference: loading this file explicitly requests the keto target.

The implementation imposes existing bond orders and formal charges; it does
not add/remove atoms or bonds, search for a compatible Kekule form outside the
matched core, or validate and roll back products. These templates include the
whole fused system to avoid leaving an adjoining aromatic ring inconsistent.
Do not apply that reasoning to a partial fused system.

Validate changes with `tstandardise` (check its diagnostics and failure file,
not just its exit status), formula/valence checks, and representative controls.
The executable regression is `src/test/fileconv/external_guanine` and requires
`fileconv` and `tstandardise`, without Python bindings or RDKit.

## Hypoxanthine

Load `hypoxanthine.textproto` after the built-in rules:

```shell
fileconv -g all -g EXT:${LILLYMOL_HOME}/data/standardisation/hypoxanthine.textproto \
  -o smi -S standardised input.smi
```

This specification chooses a whole-core keto representation for supported
hypoxanthines, including alternative keto N-H placements. It overrides the
built-in result wherever a query matches, rather than merely converting any
OH left behind. In an unsubstituted core the target is the 1H,9H keto form:

```text
OC1=NC=NC2=C1N=CN2   -> O=C1NC=NC2=C1N=CN2
O=C1NC=NC2=C1NC=N2   -> O=C1NC=NC2=C1N=CN2
O=C1N=CNC2=C1N=CN2   -> O=C1NC=NC2=C1N=CN2
```

Four mutually exclusive queries handle combinations of two choices:

- When N3 is unsubstituted, impose the N1-amide/N3-imine pattern. N1 can
  retain a substituent. When N3 is substituted, keep N1 unsubstituted and
  impose the alternative pattern in which N3 has three single bonds.
- When N7 is unsubstituted, impose the N9-pyrrole pattern (9H if N9 is
  unsubstituted). When N7 is substituted, keep N9 unsubstituted and impose
  the compatible N7-pyrrole pattern instead.

The whole matched core must be neutral and aromatic, in a fused system of
exactly two rings. C2 and C8 must each have two heavy-atom connections and
one hydrogen. O6 must be terminal. More extensively fused or reduced cores,
C2/C8 substitution, sulfur replacements and O-alkyl groups are outside this
initial scope. Simultaneous N1/N3 substitution or N7/N9 substitution is also
excluded. Isotope labels and remote stereochemistry are retained. Explicit
hydrogen atoms are not moved; use implicit H at transfer sites.

The query atom-map numbers annotate purine positions (O6 is labelled 10).
LillyMol ignores these maps for matching and for template correspondence:
the target SMILES still follows the query's positional atom order, listed
in the query-file comments. Keep all five hypoxanthine files together.

As with guanine, loading this file explicitly requests its keto/N-H target,
even with `-g all -g ENOL`. Guanine and hypoxanthine specifications can be
loaded together; their different C2 requirements keep their scopes separate.
Neither specification is enabled by `-g all` alone.

The executable regression is `src/test/fileconv/external_hypoxanthine`.
Validation on the raw ChEMBL hypoxanthine-family candidates found one
pre-existing repeat-processing failure with built-in `-g all`, CHEMBL83601.
Its C8-amino group is outside these templates, so this specification does
not resolve that failure. A matched query does not imply that every possible
representation of the molecule is covered by the current scope.
