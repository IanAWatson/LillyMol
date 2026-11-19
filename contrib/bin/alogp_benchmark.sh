#!/usr/bin/env bash
# Wrapper for alogp_benchmark.py. Needs rdkit and numpy, and the LillyMol python
# bindings, so it runs through run_python.sh. Set PYTHON to an interpreter that
# matches lib/lillymol.so.soabi and has rdkit installed.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
: "${LILLYMOL_HOME:=$(cd -- "${script_dir}/../.." && pwd)}"
export LILLYMOL_HOME
exec "${LILLYMOL_HOME}/run_python.sh" "${script_dir}/alogp_benchmark.py" "$@"
