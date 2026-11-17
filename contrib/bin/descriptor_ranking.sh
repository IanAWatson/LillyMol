#!/usr/bin/env bash
# Wrapper for descriptor_ranking.py. The script needs numpy, pandas, scipy,
# scikit-learn and matplotlib, but not the LillyMol python bindings; set PYTHON
# to an interpreter that has them.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
exec "${PYTHON:-python3}" "${script_dir}/descriptor_ranking.py" "$@"
