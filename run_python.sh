#!/usr/bin/env bash

# Run Python with PYTHONPATH and LD_LIBRARY_PATH pointing to the shared
# libraries staged by src/copy_shared_libraries.sh.

PYTHON=${PYTHON:-python3}
EXPECTED_VERSION="3.13"
actual_version=$(${PYTHON} -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
if [[ "${actual_version}" != "${EXPECTED_VERSION}" ]]; then
  echo -e "\nWarning: python version mismatch: got ${actual_version} need ${EXPECTED_VERSION}. Likely problems ahead\n" >&2
fi

my_dir=$(dirname "$0")

if [[ -v LILLYMOL_HOME ]] ; then
  BDB="${LILLYMOL_HOME}/third_party/BDB/lib"
else
  BDB="${my_dir}/third_party/BDB/lib"
fi

if [[ ! -d ${BDB} ]] ; then
  BDB=""
fi

export LD_LIBRARY_PATH=${my_dir}/lib:${BDB}:${LD_LIBRARY_PATH}
python_flags=()
if ${PYTHON} -h 2>&1 | grep -q -- " -P "; then
  python_flags=(-P)
fi

PYTHONPATH=${my_dir}/lib:${my_dir}/src:${PYTHONPATH} ${PYTHON} "${python_flags[@]}" "$@"
