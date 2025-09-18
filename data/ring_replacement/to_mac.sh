#!/bin/bash

# Rename hidden files to case insensitive forms

for file in $(ls rings_*.smi) ; do
  newname=$(echo $file | sed -e 's/A/Al/g' -e 's/a/Ar/g')
  mv $file $newname
done
