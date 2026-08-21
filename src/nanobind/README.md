# Nanobind pilot

This directory is an experimental nanobind implementation of LillyMol Python
bindings. It is intentionally built in parallel with `src/pybind` so the current
pybind11 bindings remain available while migration issues are explored.

The long-term goal is a complete changeover, not a mixed pybind11/nanobind
runtime where Python code passes C++ objects between binding systems. At
changeover, the nanobind extension will become the normal `lillymol` Python
module and the temporary `lillymol_nb` module name will disappear. During the
pilot, nanobind modules should therefore be designed as if they will eventually
own the Python-visible LillyMol types. Avoid depending on pybind-bound classes.

`Atom`, `Bond`, and `Ring` objects returned from a `Molecule` are borrowed views
into that molecule. They are valid while the parent molecule remains alive and
structurally unchanged. Structural mutation of the parent molecule can invalidate
previously returned child objects; the pilot does not add extra ownership or
copying solely to guard against that case. APIs that return independent objects
should do so explicitly, as `rings()` does by returning copied rings.

Current pilot module:

```shell
bazel build -c opt //nanobind:lillymol_nb
bazel test -c opt //nanobind:lillymol_nb_test
```

The binding implementation is split by functional area:

- `lillymol_nb.cc` contains only the `NB_MODULE` entry point;
- `lillymol_nb_internal.h` contains shared declarations for the pilot;
- `lillymol_nb_common.cc` contains shared helper implementations;
- `lillymol_nb_io.cc`, `lillymol_nb_molecule.cc`,
  `lillymol_nb_atom_bond.cc`, `lillymol_nb_set_of_atoms.cc`,
  `lillymol_nb_substructure.cc`, `lillymol_nb_descriptors.cc`,
  `lillymol_nb_standardise.cc`, and `lillymol_nb_fingerprint.cc` contain the
  Python-visible bindings.

To stage the nanobind modules separately from the pybind modules:

```shell
./copy_nanobind_shared_libraries.sh
PYTHONPATH=/path/to/LillyMol/nanobind_lib python3
```

The default destination is the repository-level `nanobind_lib` directory when run
from `src`. Pass an explicit directory to override it. Keep this directory
separate from `lib` while both pybind and nanobind builds exist.

The Python 3.13 interpreter installed by `uv` is externally managed, so install
extra packages into a virtual environment rather than modifying that base
interpreter. One working setup is:

```shell
uv venv --python /home/ian/.local/bin/python3.13 /home/ian/lillymol_py313_venv
uv pip install --python /home/ian/lillymol_py313_venv/bin/python numpy protobuf
```

`lillymol_nb` currently validates:

- constructing, copying, indexing, and iterating over a `Molecule`;
- common `Molecule` accessors, RDKit-style aliases, mutators, bond lookup,
  formal-charge, isotope, Hydrogen, ring, aromaticity, symmetry, path, SMILES,
  atom-renumbering, organic-element, and chirality helpers;
- `MolFromSmiles` returning an owned `Molecule` or `None`;
- basic string, scalar, tuple, optional, and `list[int]` conversions;
- molecule-owned `Atom` and `Bond` references from `atom(i)`, `mol[i]`,
  `for atom in mol`, `atom[i]`, `for bond in atom`, `bond(i)`, and `bonds()`;
- copied `Chiral_Centre` values, `ChiralType`, `is_actually_chiral`, and
  `tetrahedral_chirality`;
- borrowed `Atom` and `Bond` views reflecting parent molecule property mutations;
- sequence-like `Set_of_Atoms` and `Ring` atom-number iteration;
- `Ring` bindings, including both molecule-owned references from `ring(i)` and
  copied rings from `rings()`;
- `Molecule.__contains__` overloads for atomic numbers, element symbols, and
  `SubstructureQuery`;
- `SubstructureQuery` and `SubstructureResults` bindings with copied embeddings
  and borrowed-result iteration via `for embedding in sresults`;
- `TSubstructure` query-list workflows, including batch search, match counts,
  query modifiers, matched-atom reporting, and isotope labelling;
- transient SMARTS convenience helpers (`HasSubstructMatch`,
  `CountSubstructMatches`, `GetSubstructMatches`);
- GIL release around substructure search calls;
- simple descriptor helper dependencies (`alogp`, `xlogp`, `tpsa`);
- basic `Chemical_Standardisation` access via `Standardise`;
- list-returning and NumPy-returning linear, ECFP, and atom-pair fingerprint
  helpers plus `tanimoto`;
- basic molecule file I/O via `Reader`, `Writer`, `MolReaderContext`, and
  `MolWriterContext`, with `ReaderContext` and `ContextWriter` retained as
  compatibility aliases;
- `slurp` and the small global MDL/SDF reader option helpers exposed by the
  pybind I/O module.

See `changeover_checklist.md` for the current pybind-to-nanobind gap triage and
changeover tasks.

Still to explore:

- broader query functionality, including proto construction and richer query-file
  workflows;
- broader object lifetime and reference return policies for molecule-owned
  objects, especially if mutable atom or bond operations are added;
- final packaging work for replacing the pybind11 modules with nanobind modules
  under the public `lillymol` name.
