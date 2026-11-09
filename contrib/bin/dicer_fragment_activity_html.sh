#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
export LILLYMOL_HOME="${LILLYMOL_HOME:-$(cd -- "${script_dir}/../.." && pwd)}"
# Use the checked LillyMol environment when staged bindings are available.
# Plain JSON-to-HTML rendering still works without a LillyMol installation.
if [[ -f "${LILLYMOL_HOME}/lib/lillymol.so" ]]; then
  exec "${LILLYMOL_HOME}/run_python.sh" "${script_dir}/dicer_fragment_activity_html.py" "$@"
fi
exec "${PYTHON:-python3}" "${script_dir}/dicer_fragment_activity_html.py" "$@"
