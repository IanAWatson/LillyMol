#!/bin/bash

here=$(dirname $0)
fname=$(basename $0)
exec ${here}/../../bin/$(uname)/${fname%%.sh} "$@"
