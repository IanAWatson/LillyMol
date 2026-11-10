#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
export LILLYMOL_HOME="${LILLYMOL_HOME:-$(cd -- "${script_dir}/../.." && pwd)}"
if [[ -f "${LILLYMOL_HOME}/lib/lillymol.so" ]]; then
  exec "${LILLYMOL_HOME}/run_python.sh" "${script_dir}/descriptor_exploration.py" "$@"
fi
exec "${PYTHON:-python3}" "${script_dir}/descriptor_exploration.py" "$@"
