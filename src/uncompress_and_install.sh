#!/usr/bin/env bash

# Parts of LillyMol come in compressed form or with obfuscated file names.
# Convert those.
# Compressed sources are retained so installation can be repeated.

dicer_fragments_dir="${LILLYMOL_HOME}/data/dicer_fragments/"
if [[ -s "${dicer_fragments_dir}/dicer_fragments.tar.xz" ]] ; then
  echo "Uncompressing fragments"
  (cd ${dicer_fragments_dir} && tar Jxf dicer_fragments.tar.xz)
fi

chembl_sidechains="${LILLYMOL_HOME}/data/chembl_sidechains.textproto.xz"
if [[ -s "${chembl_sidechains}" ]] ; then
  echo "Uncompressing CHEMBL sidechains"
  (
    # Stage beside the destination so the final rename is atomic. Preserve
    # both the compressed source and any existing output if decompression fails.
    temporary=$(mktemp "${chembl_sidechains%.xz}.XXXXXX") || exit 1
    trap 'rm -f -- "$temporary"' EXIT
    xz --decompress --stdout -- "${chembl_sidechains}" > "$temporary" || exit 1
    mv -f -- "$temporary" "${chembl_sidechains%.xz}"
  ) || exit $?
fi

ring_replacement="${LILLYMOL_HOME}/data/ring_replacement/"

if [[ -d "${ring_replacement}" ]] ; then
  "${PYTHON:-python3}" "${ring_replacement}/install.py" || exit $?
fi
