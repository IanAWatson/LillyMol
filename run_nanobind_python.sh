#!/usr/bin/env bash

# Compatibility wrapper from the nanobind pilot. The standard hard-switched
# runtime is now run_python.sh with libraries staged in ../lib.

script_dir=$(cd "$(dirname "$0")" && pwd)
exec "${script_dir}/run_python.sh" "$@"
