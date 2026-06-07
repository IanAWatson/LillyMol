#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

html="$(mktemp /tmp/smiles_browser.XXXXXX.html)"

python3 "${SCRIPT_DIR}/smiles_browser.py" "$@" -o "${html}"

if [[ ! -s "${html}" ]] ; then
  echo "Did not create html" >&2
  exit 1
fi

xdg-open "${html}"
