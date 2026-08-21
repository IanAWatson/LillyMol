# Nanobind Changeover Checklist

This checklist tracks remaining pybind-to-nanobind compatibility work before the
nanobind bindings replace the public `lillymol` Python module name.

## Completed Baselines

- Core molecule construction, copying, indexing, iteration, atom/bond/ring views.
- RDKit-style aliases for common molecule, atom, bond, ring, chirality, and
  substructure operations.
- RingInfo-style helper object.
- SubstructureQuery, SubstructureResults, TSubstructure, and transient SMARTS
  convenience functions.
- File I/O contexts and compatibility aliases for context-aware reader/writer
  names.
- Fingerprint helpers, including list-returning and NumPy-returning paths.
- Descriptor and tool helpers: ALogP, XLogP, TPSA, QED, IWDescr, JWCats,
  MFormula, donor/acceptor assignment, charge assignment, Dicer,
  ring replacement, MedchemWizard, truncated distance matrix, GFP context, GFP
  server, reactions, and optional BerkeleyDB StructureDatabase bindings.
- Recently added Molecule helper parity: fragment membership/removal,
  largest-fragment reduction, component creation, graph/scaffold conversion,
  atom sorting, distance-matrix recomputation, atom/ring-system labelling,
  SMARTS/aromatic helpers, isotope-vector setting, directional-bond cleanup, and
  connection-table reordering.
- Python 3.13 nanobind test environment with NumPy and protobuf:
  `/home/ian/lillymol_py313_venv`.

## Remaining Pybind-only Symbols From Source Scan

The current crude `.def(...)` name comparison reports these pybind-only names.
Some are intentionally low priority; this list is for triage, not an automatic
porting mandate.

| Symbol | Area | Status |
| --- | --- | --- |
| `construct_from_proto` | Query/proto | Implemented in nanobind by accepting a Python protobuf object, calling `SerializeToString()`, and parsing the C++ proto before `ConstructFromProto`. This is appropriate for configuration/setup use; revisit only if high-volume proto transfer becomes important. |
| `rxn_from_file` | Reaction convenience | Likely easy wrapper if existing Python code imports it. Triage for compatibility. |
| `dihedral_scan` | Coordinate/conformer helper | Deferred unless users need this from Python; user explicitly does not want broad conformer-search work. |
| `non_sssr_ring` | Ring detail | Implemented in nanobind. |
| `non_sssr_rings` | Ring detail | Implemented in nanobind. |
| `fused_neighbours` | Ring detail | Possibly useful compatibility alias; check nanobind Ring already exposes enough. |
| `GetNeighbors` | Atom/RDKit alias | Low effort if expected by RDKit-oriented callers. |
| `set_vector` | Set_of_Atoms mutability | Low priority legacy/container helper. |
| `scatter` | Set_of_Atoms vector helper | Low priority unless existing scripts use it. |
| `increment_vector` | Set_of_Atoms vector helper | Low priority unless existing scripts use it. |
| `clear` | container helpers | Low priority; Python list-like behavior may be enough. |
| `pop_back` | container helpers | Low priority; Python list-like behavior may be enough. |
| `compute_Del_Re_partial_charges` | Partial charges | Disabled/commented in pybind source; do not port unless C++ support is confirmed. |
| `compute_Pullman_partial_charges` | Partial charges | Disabled/commented in pybind source; do not port unless C++ support is confirmed. |

## Changeover Work Still Needed

- Rename/package nanobind so the public module is `lillymol`; remove temporary
  `lillymol_nb` naming from installed usage.
- Decide whether to keep separate compatibility modules such as
  `lillymol_query`, `lillymol_fingerprint`, and `lillymol_tools`, or provide all
  functionality through `lillymol` with optional import aliases.
- Update installation/copy scripts so stale ABI3 artifacts cannot shadow freshly
  built nanobind modules during direct `PYTHONPATH=bazel-bin/nanobind` testing.
- Run representative user scripts against nanobind-only imports before removing
  pybind from the active path.
