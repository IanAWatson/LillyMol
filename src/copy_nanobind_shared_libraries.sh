#!/usr/bin/env bash

# Compatibility wrapper from the nanobind pilot. The production staging script is
# now copy_shared_libraries.sh, which stages nanobind modules by default.

set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo_home=$(cd "${here}/.." && pwd)

if [[ $# -gt 0 ]] ; then
  destdir=$1
else
  destdir="${repo_home}/nanobind_lib"
fi

exec "${here}/copy_shared_libraries.sh" "${destdir}"
