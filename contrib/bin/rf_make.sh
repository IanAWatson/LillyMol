#!/usr/bin/env bash

source "$(dirname -- "${BASH_SOURCE[0]}")/lib/model_environment.sh" || exit 1
check_model_protos || exit 1
exec "${PYTHON:-python}" "${model_bin_dir}/xgbd/rf_make.py" "$@"
