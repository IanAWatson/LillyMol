#!/usr/bin/env bash

# Parts of LillyMol come in compressed form or with obfuscated file names.
# Convert those.
# If this script has been run and succeeded it will be a no-op

dicer_fragments_dir="${LILLYMOL_HOME}/data/dicer_fragments/"
if [[ -s "${dicer_fragments_dir}/dicer_fragments.tar.xz" ]] ; then
  echo "Uncompressing fragments"
  cd ${dicer_fragments_dir} && tar Jxf dicer_fragments.tar.xz
fi

ring_replacement="${LILLYMOL_HOME}/data/ring_replacement/"

if [[ -d "${ring_replacement}" ]] ; then
  if [[ $(uname) == "Linux" ]] ; then
    cd ${LILLYMOL_HOME}/data/ring_replacement && ./to_linux.sh
  elif [[ $(uname) == "Darwin" ]] ; then
    cd ${LILLYMOL_HOME}/data/ring_replacement && ./to_mac.sh
  fi
fi
