#!/usr/bin/env bash

# Run python with PYTHONPATH pointing to the shared libraries generated
# by build_from_src.sh
# Make sure that PATH is set so that the same version of python that
# is in MODULE.bazel is used for run-time.
# This will not work otherwise.

# build_from_src.sh must have been run, and one of the last things it
# does is to run copy_shared_libraries.sh, which copies the compiled
# shared libraries out of bazel-bin to ../lib

EXPECTED_VERSION="3.13"
actual_version=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
if [[ "${actual_version}" != "${EXPECTED_VERSION}" ]]; then
  echo -e "\nWarning: python version mismatch: got ${actual_version} need ${EXPECTED_VERSION}. Likely problems ahead\n" >&2
fi

# Assume that lib, with the shared libraries, is at the same level as this script.
my_dir=$(dirname $0)
export LD_LIBRARY_PATH=${my_dir}/lib:${LD_LIBRARY_PATH}
PYTHONPATH=${my_dir}/lib:${my_dir}/src:${PYTHONPATH} python3 -P "$@"
