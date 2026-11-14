#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
export LILLYMOL_HOME="${LILLYMOL_HOME:-$(cd -- "${script_dir}/../.." && pwd)}"
export PYTHONPATH="${script_dir}:${PYTHONPATH:+:${PYTHONPATH}}"
if [[ -f "${LILLYMOL_HOME}/lib/lillymol.so" ]]; then
  exec "${LILLYMOL_HOME}/run_python.sh" "${script_dir}/descriptor_exploration.py" "$@"
fi
# Observe that when we run this on large multi-core machines
# we can see severely diminished performance.
# Codex suggested this fix.
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
exec "${PYTHON:-python3}" "${script_dir}/descriptor_exploration.py" "$@"
