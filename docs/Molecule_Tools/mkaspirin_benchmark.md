# Recursive aspirin SMARTS: Python binding timings

Measured 18 November 2026. Query generation is described in [mkaspirin](mkaspirin.md).

Intel Core i7-8700K (3.70 GHz), Linux, Python 3.12, RDKit 2024.09.5 and staged
LillyMol nanobind. Both engines perform first-hit searches in one thread.
Imports, file reads, molecule construction/perception and query generation are
outside the measured sections. Each row is the median of three fresh Python
worker processes. Warm-search measurements are medians of five subsequent
searches per process. Panel measurements are medians of three 11-molecule passes
per process. All expected matches were checked on every panel pass.

The same portable SMARTS bytes and molecule SMILES were supplied to both engines.
Each query is generated with seed 19, one recursive alternative, default contexts,
maximum depth 32. Byte counts exclude names and newlines. The panel contains
aspirin, aspirin with a disconnected sodium ion, and nine related negatives.

All times below are milliseconds. RSS is MiB and includes the Python process,
libraries, query, molecules and matching allocations, not just the parsed query.

| Bytes | Depth | LillyMol parse | RDKit parse | LillyMol first search | RDKit first search | LillyMol warm search | RDKit warm search | LillyMol 11-molecule panel | RDKit panel | LillyMol peak RSS | RDKit peak RSS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1144 | 2 | 0.892 | 0.454 | 0.168 | 0.106 | 0.049 | 0.019 | 0.273 | 0.155 | 57.9 | 60.1 |
| 10061 | 6 | 13.830 | 12.263 | 2.347 | 0.799 | 0.780 | 0.471 | 2.428 | 3.536 | 57.9 | 61.6 |
| 100086 | 11 | 42.326 | 104.274 | 10.511 | 5.266 | 5.095 | 5.028 | 10.388 | 49.259 | 64.1 | 76.3 |
| 1000141 | 15 | 628.326 | 1991.523 | 100.234 | 63.047 | 51.194 | 60.181 | 97.530 | 571.099 | 430.2 | 227.3 |
| 3000360 | 19 | 1781.517 | 7244.209 | 295.135 | 196.970 | 141.260 | 173.741 | 283.024 | 1695.810 | 1244.1 | 559.3 |

## LillyMol-extension queries

These use different random expansion sequences and predicate lengths from the
portable queries. Equal file sizes are not a controlled measure of extension
overhead. RDKit was not asked to parse LillyMol extensions.

| Bytes | Depth | Parse ms | First search ms | Warm search ms | Panel ms | Peak RSS MiB |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1098 | 2 | 0.321 | 0.075 | 0.021 | 0.120 | 57.9 |
| 10016 | 8 | 17.639 | 1.982 | 0.935 | 2.420 | 57.9 |
| 100492 | 16 | 137.724 | 8.624 | 4.311 | 7.456 | 59.8 |
| 1000510 | 22 | 521.492 | 77.838 | 41.507 | 83.590 | 382.9 |
| 3000029 | 25 | 1828.716 | 233.656 | 123.238 | 247.099 | 1101.2 |

These are deliberately pathological aspirin queries and a small test panel,
not a general comparison of toolkit substructure performance. Results depend
on query shape, target molecules, versions, build options and hardware.

## Reproduce the measurements

The [benchmark driver](benchmarks/recursive_aspirin.py) generates both query
families, runs each configuration in three fresh processes, verifies the panel
matches and writes process-level measurements to `results.json`. Query files
are generated locally rather than checked into the repository.

Run through LillyMol's wrapper using an environment with matching staged
bindings and RDKit installed:

```shell
PYTHON=/path/to/venv/bin/python3 \
  ${LILLYMOL_HOME}/run_python.sh \
  ${LILLYMOL_HOME}/docs/Molecule_Tools/benchmarks/recursive_aspirin.py \
  --output_dir /tmp/aspirin_timings
```

Each worker has a 60-second timeout. Errors and timeouts are recorded rather
than treated as timing results. The existing output files in the specified
directory are replaced, so use a dedicated directory. RSS units are normalized
to MiB on Linux and macOS.

The recorded [raw measurements](benchmarks/recursive_aspirin_2026-11-18.json)
include every process-level result. Stored query paths are reduced to filenames;
the driver reproduces those files in its output directory.
