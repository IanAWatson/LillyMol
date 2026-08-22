#!/usr/bin/env bash

# Run python with PYTHONPATH and LD_LIBRARY_PATH pointing to the nanobind
# shared libraries staged by src/copy_nanobind_shared_libraries.sh.
# Make sure that PATH is set so that the same version of python that
# is in MODULE.bazel is used for run-time.
# This will not work otherwise.

# Build the nanobind targets first, then run src/copy_nanobind_shared_libraries.sh
# to copy extension modules and dependent shared libraries to ../nanobind_lib.

PYTHON=${PYTHON:-python3}

EXPECTED_VERSION="3.13"
actual_version=$(${PYTHON} -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
if [[ "${actual_version}" != "${EXPECTED_VERSION}" ]]; then
  echo -e "\nWarning: python version mismatch: got ${actual_version} need ${EXPECTED_VERSION}. Likely problems ahead\n" >&2
fi

# Assume that nanobind_lib, with the shared libraries, is at the same level as this script.
my_dir=$(dirname $0)

# BerkekelyDB libs, may or may not be present, was BUILD_BDB set during the build.
if [[ -v LILLYMOL_HOME ]] ; then
  BDB="${LILLYMOL_HOME}/third_party/BDB/lib"
else
  BDB="${my_dir}/third_party/BDB/lib"
fi

# Perhaps BDB not installed.
if [[ ! -d ${BDB} ]] ; then
  BDB=""
fi

export LD_LIBRARY_PATH=${my_dir}/nanobind_lib:${BDB}:${LD_LIBRARY_PATH}
python_flags=()
if ${PYTHON} -h 2>&1 | grep -q -- " -P "; then
  python_flags=(-P)
fi

PYTHONPATH=${my_dir}/nanobind_lib:${my_dir}/src:${PYTHONPATH} ${PYTHON} "${python_flags[@]}" "$@"
