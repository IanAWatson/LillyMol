# Built-in LillyMol Binding Benchmark Results

These are small smoke-benchmark results for the pybind11 and nanobind Python
bindings. They use the 10 built-in molecules in
`benchmark_lillymol_bindings.py`, so the numbers are useful for quick local
comparison rather than final performance claims.

## Run Context

- Date: 2026-08-21
- Source checkout HEAD at run time: `dc7aa79d`
- Python: `3.13.7`, `/home/ian/.local/bin/python3.13`
- Platform: `Linux-6.8.0-137-generic-x86_64-with-glibc2.39`
- Molecules: 10 built-in SMILES
- Loops per benchmark: 20
- Raw pybind results: `results_builtin_pybind.json`
- Raw nanobind results: `results_builtin_nanobind.json`

Commands:

```shell
PYTHONPATH="$PWD/bazel-bin/pybind:$PWD/bazel-bin" \
LILLYMOL_HOME=/home/ian/LillyMol_IWFORK \
/home/ian/.local/bin/python3.13 \
  nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding pybind --loops 20 --json \
  > nanobind/benchmarks/results_builtin_pybind.json

PYTHONPATH="$PWD/bazel-bin/nanobind:$PWD/bazel-bin" \
LILLYMOL_HOME=/home/ian/LillyMol_IWFORK \
/home/ian/.local/bin/python3.13 \
  nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding nanobind --loops 20 --json \
  > nanobind/benchmarks/results_builtin_nanobind.json
```

## Common Benchmarks

Lower `best_s` is faster. The ratio is `nanobind best_s / pybind best_s`.

| Benchmark | pybind best s | nanobind best s | Ratio |
| --- | ---: | ---: | ---: |
| atom_index_access | 0.000063839 | 0.000023700 | 0.37 |
| atom_iteration | 0.000180625 | 0.000053225 | 0.29 |
| parse | 0.000297113 | 0.000094976 | 0.32 |
| scalar_properties | 0.000010599 | 0.000002525 | 0.24 |
| substructure_count | 0.000107317 | 0.000099984 | 0.93 |
| unique_smiles | 0.000004447 | 0.000002784 | 0.63 |

## Skips

NumPy was not installed in this Python 3.13 environment. NumPy-specific
benchmarks were therefore skipped. In addition, importing the pybind
`lillymol_fingerprint` module in this environment required NumPy, so pybind
fingerprint benchmarks were skipped in this run.

Nanobind-only benchmarks that completed in this run include
`substructure_embeddings`, `linear_fingerprint`, `ecfp`, and `coordinates_list`;
see `results_builtin_nanobind.json` for timings.
