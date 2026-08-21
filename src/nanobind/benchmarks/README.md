# LillyMol Python Binding Benchmarks

These benchmarks compare the pybind11 and nanobind LillyMol Python bindings on
small Python-facing workflows. They are timing tools, not pass/fail regression
tests.

Build the bindings first from `src`:

```shell
bazel build -c opt //pybind:lillymol //pybind:lillymol_fingerprint //pybind:lillymol_tools //pybind:lillymol_query
bazel build -c opt //nanobind:lillymol_nb
```

Run pybind:

```shell
PYTHONPATH="$(pwd)/bazel-bin/pybind:$(pwd)/bazel-bin" \
  python3 nanobind/benchmarks/benchmark_lillymol_bindings.py --binding pybind
```

Run nanobind:

```shell
PYTHONPATH="$(pwd)/bazel-bin/nanobind:$(pwd)/bazel-bin" \
  python3 nanobind/benchmarks/benchmark_lillymol_bindings.py --binding nanobind
```

For a useful ChEMBL-derived input set, a few thousand molecules is enough for
quick iteration and 50k-100k molecules is useful for steadier whole-workflow
measurements. Use a normal LillyMol SMILES file, one molecule per line:

```shell
python3 nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding nanobind \
  --smiles-file /path/to/random_chembl.smi \
  --max-mols 10000
```

Descriptor benchmarks require `LILLYMOL_HOME` to point at a tree containing the
standard query data. NumPy-specific benchmarks are skipped if NumPy is not
installed or if the selected binding does not expose the needed method.

Checked-in baseline runs:

- `results_builtin.md`: 10 built-in molecules, useful as a smoke benchmark.
- `results_rand10k.md`: 10,000 random molecules from `/home/ian/rand10k.smi`,
  useful as the current larger comparison baseline.

Raw JSON is stored beside each summary.
