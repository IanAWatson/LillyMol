#!/bin/bash
here=$(dirname $0)
if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(dirname $here))
fi
export PATH=${here}:$PATH
exec perl ${here}/gfp_make.pl -bindir ${here} "$@"
