# Agent Notes for LillyMol

This file gives AI/code agents a compact orientation to the LillyMol repository.
It is not a replacement for reading the relevant source, tests, and docs. Treat
it as starting context and verify details locally before changing behavior.

## Repository layout

- `src/` contains most C++ source and is usually the working directory.
- `src/Molecule_Lib/` contains the core molecule representation, SMILES,
  substructure searching, atom typing, and low-level chemistry functionality.
- `src/Molecule_Tools/` contains many command-line molecular tools.
- `src/Utilities/GFP_Tools/` contains fingerprint and similarity-search tools.
- `src/Molecule_Tools_Bdb/` contains BerkeleyDB-backed tools.
- `src/pybind/` contains Python bindings built with pybind11.
- `docs/` contains user documentation. Tool docs are commonly under
  `docs/Molecule_Tools`, `docs/GFP`, `docs/Molecule_Lib`, and `docs/python`.
- `test/` contains executable regression tests driven by `test/run_all_test.rb`.
- `data/` contains query files, standardization resources, and other runtime
  data. Many tools expect `LILLYMOL_HOME` to point at the repository root or an
  installed LillyMol tree.

## Build and test conventions

- Prefer Bazel for builds. Example from `src/`:

  ```shell
  bazel build //Molecule_Tools:common_names
  ```

- Executable tests are usually run with:

  ```shell
  ruby ../test/run_all_test.rb -bindir /path/to/bin tool_name
  ```

  When using Bazel-built tools from `src/`, a typical `-bindir` is
  `/path/to/LillyMol/src/bazel-bin/Molecule_Tools`.

- Most executable tests live under `test/<tool>/<case>/` with:
  - `in/` for inputs;
  - `out/` for expected outputs;
  - `tests.json` at the tool directory level.

- Tests commonly compare files produced by `-S`, `-T`, or similar output-stem
  options. If a file exists in the expected `out/` directory, the test runner can
  compare the generated file of the same name.

- Some tests use a shell script in the test case directory when the standard
  command-line construction is insufficient. Keep scripts deterministic and make
  sure they run in the temporary test directory supplied by the runner.

- Add focused regression tests when changing behavior. LillyMol has many older
  tools with broad option surfaces; small tests are usually more valuable than
  large end-to-end fixtures.

## General coding conventions

- Follow the style of nearby code. Much of LillyMol predates modern C++; avoid
  broad modernization inside unrelated changes.
- Prefer existing LillyMol types and utilities in older code paths:
  `IWString`, `const_IWSubstring`, `IWString_and_File_Descriptor`,
  `resizable_array`, `resizable_array_p`, `Set_of_Atoms`, and related classes.
- Use `rg` for source searches.
- Preserve command-line compatibility unless the requested change explicitly
  allows breaking behavior.
- Be careful with file-scope static/global state. Some older tools use it
  heavily; avoid introducing new global mutable state when an options/class
  object would be straightforward.
- Many tools support `-` as stdin/stdout. Preserve that convention where it
  already exists.
- Many molecule-writing tools use `Molecule_Output_Object`. Before opening
  output from a user-supplied stem, check whether it would overwrite input files,
  typically with `would_overwrite_input_files`.
- Two-pass tools are especially vulnerable to input clobbering: opening an
  output file before the second input pass can truncate the file that still needs
  to be read.
- `IWString_and_File_Descriptor` is widely used for buffered text output. A file
  name of `-` often means stdout, but confirm the specific code path.
- Avoid destructive git operations. The worktree often contains user changes and
  generated files.

## Documentation conventions

- User documentation should start with what the tool does and how to use it.
  Put build/setup/background material later unless it is required for first use.
- Prefer concrete command examples with short input/output snippets.
- Document risk/behavior changes explicitly when an option broadens chemical
  equivalence, changes output semantics, or can generate many more structures.
- Keep option documentation tied to observed code behavior. Some older docs may
  describe stale or partially implemented options.
- If a documented command is practical to run as a regression test, consider
  adding it under `test/<tool>`.

## Chemistry and LillyMol behavior notes

- LillyMol unique SMILES is often used as the structural identity key, but many
  tools deliberately transform molecules before generating that key.
- Representative-output tools often write the first input structure encountered
  for a group, even when a transformed molecule was used for comparison.
- Standardization, largest-fragment selection, isotope removal, graph-form
  comparison, and element transformations can all change equivalence classes.
- Atom typing is compositional. Bit masks in `Molecule_Lib/atom_typing.h` can be
  combined into a `uint32_t` atom type.
- SMARTS/substructure code includes LillyMol-specific extensions such as global
  ids (`/IWgid...`) and proto-based query features. Verify behavior in
  `Molecule_Lib` before documenting or changing it.
- Ring and ring-system substituent matching is greedy by design; ordering of
  query specifications can matter.
- GFP/TDT tools commonly use records such as `$SMI<`, `PCN<`, and `DIST<`, with
  items separated by `|`.
- Many pybind tests require a Python environment with the expected packages and
  shared libraries. The repository often uses wrapper scripts or `LILLYMOL_HOME`
  to locate runtime data and libraries.

## Practical workflow for agents

1. Inspect the relevant source, tests, and existing docs before proposing edits.
2. Identify whether the request is a behavior change, refactor, documentation
   update, or test addition; avoid expanding scope unnecessarily.
3. Make the smallest coherent change that preserves existing behavior outside
   the requested area.
4. Build the directly affected target when practical.
5. Run focused tests for the affected tool or library.
6. Report exactly what changed, what was tested, and any remaining caveats.
