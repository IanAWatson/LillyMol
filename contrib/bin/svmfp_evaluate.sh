#!/bin/bash
here=$(dirname $0)
if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(dirname ${here}))
fi
export PATH=${here}:${here}/../../bin/Linux:$PATH
exec ruby ${here}/svmfp_evaluate.rb "$@"
