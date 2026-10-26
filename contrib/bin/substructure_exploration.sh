#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
export LILLYMOL_HOME="${LILLYMOL_HOME:-$(cd -- "${script_dir}/../.." && pwd)}"
export PYTHONPATH="${script_dir}${PYTHONPATH:+:${PYTHONPATH}}"
exec "${LILLYMOL_HOME}/run_python.sh" "${script_dir}/substructure_exploration.py" "$@"
