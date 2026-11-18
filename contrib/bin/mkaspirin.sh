#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
export LILLYMOL_HOME="${LILLYMOL_HOME:-$(cd -- "${script_dir}/../.." && pwd)}"
exec "${LILLYMOL_HOME}/run_python.sh" "${script_dir}/mkaspirin.py" "$@"
