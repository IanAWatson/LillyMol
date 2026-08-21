# LillyMol Python wheel staging

This directory contains the minimal infrastructure needed to build a Python
wheel from Bazel-built LillyMol pybind modules.

The normal LillyMol runtime workflow is unchanged:

```shell
cd ${LILLYMOL_HOME}/src
bazel build //pybind:all
./copy_shared_libraries.sh ${LILLYMOL_HOME}/lib
```
Or if building via make,
```shell
cd ${LILLYMOL_HOME}
make all
cd src
./copy_shared_libraries.sh ${LILLYMOL_HOME}/lib
```

`run_python.sh` continues to use `${LILLYMOL_HOME}/lib` via `PYTHONPATH` and
`LD_LIBRARY_PATH`.

## Build a wheel

From `${LILLYMOL_HOME}`:

```shell
cd src
bazel build //pybind:all
./copy_shared_libraries.sh ../lib

cd ../python
./scripts/stage_wheel_files.sh
python -m pip install build wheel setuptools
python -m build --wheel
```
If using uv try
```
uv build --wheel
```
instead.

The wheel will be written to `python/dist/`.

## What gets staged

`scripts/stage_wheel_files.sh` creates `python/prebuilt/` containing:

- top-level pybind extension modules such as `lillymol.so`,
  `lillymol_query.so`, and `lillymol_tools.so`;
- compatibility extension modules such as `lillymol_atom.so` and
  `lillymol_ring.so`, which re-export classes from the canonical `lillymol`
  module;
- private LillyMol shared libraries from `${LILLYMOL_HOME}/lib`;
- selected generated protobuf Python modules and GFP HTTP helper modules.

`python/prebuilt/`, `python/build/`, `python/dist/`, and egg-info directories are
generated artifacts and are ignored by git.

## Runtime dependencies

The base wheel declares runtime dependencies on `numpy` and `protobuf`.
The GFP HTTP server helpers additionally need:

```shell
pip install 'lillymol[http]'
```

which installs `fastapi` and `uvicorn`.

## Notes

If a wheel installs but importing a module fails with a missing private shared
library, inspect the installed extension with `ldd`. Linux wheels may need a
repair step such as `auditwheel repair`, or staged extension modules may need an
`$ORIGIN` runtime path so they can find private shared libraries shipped in the
same wheel. The staging layout deliberately keeps private libraries next to the
extension modules to make that repair step straightforward.

## Detailed explanation.
I asked Codex for a more detailed information on the wheel building process,
and how it mitigated some problems. Here is the reply.

For those familiar with python wheel construction this may be non-news.



• What happened is that we turned the Bazel-built pybind outputs into a normal Python wheel and fixed the two main deployment hazards.

  First, the duplicate pybind registration problem.

  Previously, more than one extension module registered the same C++ type. For example:

  - lillymol.so registered Atom
  - lillymol_atom.so also registered Atom

  Pybind11 keeps a process-wide registry of C++ types. If Python imports two modules that both try to register the same C++ class, pybind11 can
  fail with a duplicate-registration error.

  We fixed that by making lillymol the canonical owner of core types:

  - Atom
  - Bond
  - Set_of_Atoms
  - Ring

  Then the smaller modules became shims. For example, lillymol_atom now just imports lillymol and exposes lillymol.Atom. It no longer registers
  Atom itself.

  So this is now true:

  import lillymol
  import lillymol_atom

  assert lillymol.Atom is lillymol_atom.Atom

  That removes the duplicate type registration problem while keeping old import names usable.

  Second, the shared-library loading problem.

  The pybind modules are .so files. Some are self-contained, but others depend on private LillyMol shared libraries such as:

  libiwmolecule.so
  libiwreaction.so
  libmoleculeio.so
  libdb_cxx-18.1.so
  ...

  When using run_python.sh, this works because it sets:

  LD_LIBRARY_PATH=${LILLYMOL_HOME}/lib
  PYTHONPATH=${LILLYMOL_HOME}/lib:${LILLYMOL_HOME}/src

  But a wheel cannot rely on ${LILLYMOL_HOME}/lib. It should work after pip install.

  So the wheel staging script copies:

  - top-level pybind modules like lillymol.so, lillymol_tools.so
  - private LillyMol shared libraries like libiwmolecule.so
  - generated protobuf Python files
  - selected helper Python files

  Then it uses:

  patchelf --set-rpath '$ORIGIN'

  on the staged .so files.

  $ORIGIN means “look in the same directory as this shared object”. So after the wheel is installed, lillymol_tools.so can find libdb_cxx-18.1.so
  sitting beside it inside site-packages.

  That is why the smoke test worked even with:

  env -u LD_LIBRARY_PATH

  In other words, the wheel no longer depends on run_python.sh.

  Third, the wheel packaging itself.

  The python/ directory is now a packaging wrapper around prebuilt Bazel artifacts.

  The flow is:

  cd ${LILLYMOL_HOME}/src
  bazel build //pybind:all
  ./copy_shared_libraries.sh ../lib

  cd ../python
  ./scripts/stage_wheel_files.sh
  python -m build --wheel

  stage_wheel_files.sh creates:

  python/prebuilt/

  Then setup.py copies those prebuilt files into the wheel.

  This is not compiling C++ through setuptools. Bazel still does the real build. Setuptools is just packaging already-built binaries.

  So the important conceptual split is:

  Bazel        -> builds LillyMol C++/pybind .so files
  copy script  -> places runtime .so files in LILLYMOL_HOME/lib
  stage script -> copies those into python/prebuilt and patches rpath
  setuptools   -> packages python/prebuilt into a wheel
  pip install  -> installs the wheel into site-packages

  The smoke test verified all three things:

  - Python can import every extension module from the wheel.
  - Duplicate pybind registrations are gone.
  - Private shared libraries are found without LD_LIBRARY_PATH.

  That is a solid first working wheel.

