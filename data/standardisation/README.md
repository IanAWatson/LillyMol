# Optional external chemical standardisations

These rules are opt-in additions to LillyMol chemical standardisation. They are
not enabled by `-g all` alone. For the initial guanine rules:

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
