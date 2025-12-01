#!/usr/bin/env bash

# Run all the python unit tests found in the pybind directory.

here=$(dirname $0)

if [[ ! -v PYTHONPATH ]] ; then
  export PYTHONPATH=${here}
fi

# If tmpdir is not set, multiple people running absl tests will collide.
# Unique tmpdir for each user.
if [[ ! -v TMPDIR ]] ; then
  export "TMPDIR=/tmp/absl_testing_${USER}"
fi

libdir="${here}/../lib"
if [[ ! -s "${libdir}" ]] ; then
  echo "No shared libraries available ${here}, python unit tests not done"
  exit 1
fi

if [[ ! -s "${libdir}/lillymol.so" ]] ; then
  echo "No lillymol Module ${libdir}/lillymol.so, python unit tests not done"
  exit 1
fi

if [[ ! -d "${here}/pybind" ]] ; then
  echo "pybind not found ${here}"
  exit 1
fi

run_python="${here}/../run_python.sh"
if [[ ! -x ${run_python} ]] ; then
  echo "Where is ${run_python}"
  exit 1
fi

declare -i failures=0
for file in ${here}/pybind/*_test.py ; do
  ${run_python} ${file}
  if [[ $? -ne 0 ]] ; then
    echo "${file} failed"
    failures=$(($failures + 1))
  fi
done

echo 'Python unit tests complete'
if [[ ${failures} -gt 0 ]] ; then
  echo "${failures} python tests failed"
fi
