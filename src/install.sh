#!/bin/bash

# Use bazel query to identify all cc_binary targets, and for
# those that are newer than what is in the target directory, copy
# the newer version.
# If the first argument is an existing directory, that is assumed to
# be the destination. All other arguments are passed to the `bazel query`
# command.

files_not_present=0
files_copied=0
already_copied=0

if [[ -z "$1" ]] ; then
  destination='../bin/Linux'
  echo "No argument, using default destinaton ${destination}" >&2
elif [[ -d "$1" ]] ; then
  destination="$1"
  shift
fi

# note that we can also query via tags with something like
# bazel query "attr(tags, 'tbb', //Molecule_Tools:all) intersect kind(cc_binary, //Molecule_Tools:all)" 
# But that is not available from the command line

# If this ever gets done at the directory level...
#molecule_tools=$(bazelisk query "$@" "kind(cc_binary, //Molecule_Tools:all)")
#molecule_tools_bdb=$(bazelisk query "$@" "kind(cc_binary, //Molecule_Tools_Bdb:all)")
#berkeleydb=$(bazelisk query "$@" "kind(cc_binary, //BerkeleyDB:all)")
#gfp=$(bazelisk query "$@" "kind(cc_binary, //Utilities/GFP_Tools:all)")
#gfp_knn=$(bazelisk query "$@" "kind(cc_binary, //Utilities/GFP_Knn:all)")
#general=$(bazelisk query "$@" "kind(cc_binary, //Utilities/General:all)")
#gene_expression=$(bazelisk query "$@" "kind(cc_binary, //Utilities/GeneExpression:all)")
#distance_matrix=$(bazelisk query "$@" "kind(cc_binary, //Utilities/Distance_Matrix:all)")
#tdt=$(bazelisk query "$@" "kind(cc_binary, //Foundational/iw_tdt:all)")
#obsolete=$(bazelisk query "$@" "kind(cc_binary, //Obsolete:all)")

targets=$(bazelisk query "$@" "kind(cc_binary, //...:all)")

for target in $targets ; do
  exe=$(echo ${target} | sed -e 's/..*://')
  target=$(echo ${target} | sed -e 's/\/\///' -e 's/:/\//')
  nv=bazel-bin/${target}
  ov=${destination}/${exe}
  if [[ ! -s ${nv} ]] ; then
    echo "${nv} not present"
    ((files_not_present++))
  elif [[ ${nv} -nt ${ov} ]] ; then
    echo "Copy ${exe}"
    cp -f ${nv} ${destination}
    ((files_copied++))
  else
    ((already_copied++))
  fi
done

echo "${files_not_present} files not present"
echo "${files_copied} files copied"
echo "${already_copied} files up to date"
