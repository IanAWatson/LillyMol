#!/usr/bin/env bash

export UNAME=$(uname)
export OSTYPE=$(uname)

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(cd -- "${script_dir}/../.." && pwd -P)
fi

echo "LILLYMOL_HOME ${LILLYMOL_HOME}"
export PATH=${LILLYMOL_HOME}/bin/${OSTYPE}:$PATH

exec ruby "${script_dir}/run_all_test.rb" "$@"
