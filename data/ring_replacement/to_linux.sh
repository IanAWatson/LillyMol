#!/bin/bash

# Rename hidden files to case differentiated forms

for file in $(ls rings*.smi) ; do
  newname=$(echo $file | sed -e 's/Ar/a/g' -e 's/Al/A/g')
  mv $file $newname
done
