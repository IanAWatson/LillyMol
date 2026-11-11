#!/bin/bash
# Query directories and globs given as relative paths, and a missing M: file
# alongside a good -q. $1 is the executable, $2 the 'in' directory. Runs in the
# temporary test directory, so copy the queries there and use relative paths.

exe=$1
indir=$2

cp -r "${indir}/queries" queries

echo "DIR relative"
"${exe}" -q DIR:queries -m dir_matches -o smi -i ICTE "${indir}/in.smi" 2> /dev/null
cat dir_matches.smi

echo "GLOB relative"
"${exe}" -q 'GLOB:queries/.*\.textproto' -m glob_matches -o smi -i ICTE "${indir}/in.smi" 2> /dev/null
cat glob_matches.smi

# Only nitro.textproto matches, so only nitrobenzene can match a query.
echo "GLOB selecting one file"
"${exe}" -q 'GLOB:queries/nitro.*' -m glob_one_matches -o smi -i ICTE "${indir}/in.smi" 2> /dev/null
cat glob_one_matches.smi

echo "GLOB in current directory"
cd queries
"${exe}" -q 'GLOB:.*\.textproto' -m ../glob_cwd_matches -o smi -i ICTE "${indir}/in.smi" 2> /dev/null
cd ..
cat glob_cwd_matches.smi

# A good query followed by a missing file of molecules must fail, not be ignored.
echo "missing M: file after a good query"
if "${exe}" -q queries/benzene.textproto -q M:no_such_file.smi -m ignored -i ICTE "${indir}/in.smi" 2> /dev/null; then
  echo "unexpected success"
else
  echo "failed as expected"
fi
