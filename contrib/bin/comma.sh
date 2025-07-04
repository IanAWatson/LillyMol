#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME="$(realpath $0)/../../.."
fi

QUERIES="${LILLYMOL_HOME}/data/queries"

if [[ ! -d ${QUERIES} ]] ; then
  echo "WHere are my queries ${queries}" >&2
  exit 1
fi

${LILLYMOL_HOME}/bin/Linux/comma -N noremove -N F:${QUERIES}/charges/queries \
        -H noremove -H a=F:${QUERIES}/hbonds/acceptor \
        -H d=${QUERIES}/hbonds/donor.qry -H label \
        -F ${LILLYMOL_HOME}/data/wildman_crippen.dat \
        -M min=-200 -M max=200 -b -V molvol"$@"
