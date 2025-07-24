#!/usr/bin/env bash

here=$(realpath $(dirname $0));

exec ruby ${here}/Lilly_Medchem_Rules.rb "$@"
