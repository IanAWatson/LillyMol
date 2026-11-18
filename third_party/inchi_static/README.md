# Static InChI for LillyMol

LillyMol links InChI statically, so InChI-enabled executables do not need
`libinchi.so` or `libinchi.dylib` at runtime. The normal optional build is still
`BUILD_INCHI=1 ./build_linux.sh` from `src`.

To build the dependency separately, from the repository root:

```shell
# Only needed if the checkout is absent:
git clone https://github.com/IUPAC-InChI/InChI third_party/InChI
cmake -S third_party/inchi_static -B third_party/InChI/build-static \
  -DCMAKE_BUILD_TYPE=Release
cmake --build third_party/InChI/build-static --target inchi_static --parallel 8
cd src
bazel build -c opt --config=inchi //Molecule_Tools:fileconv \
  //Molecule_Tools:unique_molecules //Molecule_Tools:common_names
```

The archive is `third_party/InChI/build-static/lib/libinchi.a`. It contains
position-independent code, so it can also be linked into shared libraries,
including Python extensions. Bazel supplies the math-library dependency.

Upstream's CMake explicitly creates a shared-library target. This wrapper loads
that target with `EXCLUDE_FROM_ALL` and reuses its source list, include paths,
compile definitions and compiler-settings dependency to create a separate static
target. It does not modify the upstream checkout or build the shared target.

The build script configures and incrementally builds the archive even when the
checkout already exists. Its separate build directory permits upgrading from
the previous shared build. No InChI library is copied into LillyMol's `lib`
directory. Old shared libraries may remain there, but newly built executables
do not depend on them; older executables can still require them.
