#!/usr/bin/env bash

source "$(dirname -- "${BASH_SOURCE[0]}")/lib/model_environment.sh" || exit 1
exec ruby "${model_bin_dir}/xgbd/xgbd_evaluate.rb" "$@"
