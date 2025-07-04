#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME="$(realpath $0)/../.."
fi

exec ${LILLYMOL_HOME}/bin/Linux/tshadow -L -G "$@"
