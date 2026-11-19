#!/usr/bin/env bash
# Built-in standardisation first; distributed external targets apply last.

if [[ $# -eq 1 && ( $1 == --help || $1 == -h ) ]]; then
  cat <<'HELP'
Usage: standardise_all.sh [fileconv options] input...

Runs fileconv with -g all and the distributed guanine and hypoxanthine
EXT specifications. External rules determine the final representation
where they match, subject to their documented substitution restrictions.

All other arguments are passed unchanged to fileconv. For example:
  standardise_all.sh -o smi -S standardised input.smi
  standardise_all.sh -l -o usmi -S - input.smi > standardised.smi
HELP
  exit 0
fi

if [[ -z ${LILLYMOL_HOME:-} ]]; then
  script_path=$(realpath "${BASH_SOURCE[0]}") || exit 1
  LILLYMOL_HOME=$(cd "$(dirname "${script_path}")/../.." && pwd) || exit 1
  export LILLYMOL_HOME
fi

specifications="${LILLYMOL_HOME}/data/standardisation"
# List top-level specifications explicitly: query files are also textprotos.
exec "${LILLYMOL_HOME}/bin/$(uname)/fileconv" \
  -g all \
  -g "EXT:${specifications}/guanine.textproto" \
  -g "EXT:${specifications}/hypoxanthine.textproto" \
  "$@"
