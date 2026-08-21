# rand10k LillyMol Binding Benchmark Results

These results use `/home/ian/rand10k.smi`, a 10,000 molecule random sample. This
is a more useful Python binding benchmark than the built-in smoke set because it
runs enough molecules for per-call overhead to show up consistently.

## Run Context

- Date: 2026-08-21
- Source checkout HEAD at run time: `dc7aa79d`
- Python: `3.13.7`, `/home/ian/.local/bin/python3.13`
- Platform: `Linux-6.8.0-137-generic-x86_64-with-glibc2.39`
- Input: `/home/ian/rand10k.smi`
- Molecules: 10,000
- Loops per benchmark: 5
- Raw pybind results: `results_rand10k_pybind.json`
- Raw nanobind results: `results_rand10k_nanobind.json`

Commands:

```shell
PYTHONPATH="$PWD/bazel-bin/pybind:$PWD/bazel-bin" \
LILLYMOL_HOME=/home/ian/LillyMol_IWFORK \
/home/ian/.local/bin/python3.13 \
  nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding pybind --smiles-file /home/ian/rand10k.smi \
  --loops 5 --json \
  > nanobind/benchmarks/results_rand10k_pybind.json

PYTHONPATH="$PWD/bazel-bin/nanobind:$PWD/bazel-bin" \
LILLYMOL_HOME=/home/ian/LillyMol_IWFORK \
/home/ian/.local/bin/python3.13 \
  nanobind/benchmarks/benchmark_lillymol_bindings.py \
  --binding nanobind --smiles-file /home/ian/rand10k.smi \
  --loops 5 --json \
  > nanobind/benchmarks/results_rand10k_nanobind.json
```

## Common Benchmarks

Lower `best_s` is faster. The ratio is `nanobind best_s / pybind best_s`.

| Benchmark | pybind best s | nanobind best s | Ratio |
| --- | ---: | ---: | ---: |
| atom_index_access | 0.168106 | 0.059509 | 0.35 |
| atom_iteration | 0.435647 | 0.101374 | 0.23 |
| parse | 0.295018 | 0.045976 | 0.16 |
| scalar_properties | 0.009007 | 0.002887 | 0.32 |
| substructure_count | 0.173931 | 0.168235 | 0.97 |
| unique_smiles | 0.004577 | 0.002363 | 0.52 |

## Nanobind-only Completed Benchmarks

These benchmarks completed for nanobind in this environment, but the comparable
pybind benchmark was unavailable or skipped.

| Benchmark | nanobind best s | Checksum |
| --- | ---: | ---: |
| coordinates_list | 0.064923 | 887082 |
| ecfp | 0.226203 | 10240000 |
| linear_fingerprint | 1.010351 | 20480000 |
| substructure_embeddings | 0.097961 | 219365 |

## Skips

NumPy was not installed in this Python 3.13 environment. NumPy-specific
benchmarks were therefore skipped. Importing the pybind `lillymol_fingerprint`
module in this environment also required NumPy, so pybind fingerprint benchmarks
were skipped in this run.
