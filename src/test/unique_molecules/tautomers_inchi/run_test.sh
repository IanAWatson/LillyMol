#! /bin/bash
# unique_molecules -C compares structures by InChI, which only exists in builds
# made with bazel --config=inchi. Check the right behaviour for whichever kind of
# build this is. Both outcomes print PASS, so a default build is not "skipped".
#
# Arguments: executable, in directory, out directory.
exe=$1
indir=$2

"$exe" -C -S unique -D duplicate -l "$indir/tautomers.smi" > c_stdout 2> c_stderr
rc=$?

if grep -q 'does not include InChI' c_stderr ; then
  # No InChI in this build: -C must be refused rather than silently ignored.
  if [ $rc -eq 0 ] ; then
    echo "FAIL: -C accepted by a build without InChI"
    exit 1
  fi
  echo "PASS"
  exit 0
fi

# Any other failure (including a missing libinchi.so) is a real failure.
if [ $rc -ne 0 ] ; then
  echo "FAIL: unique_molecules -C exited $rc"
  cat c_stderr
  exit 1
fi

# With InChI the two tautomers of 2-pyridone are the same structure, as are the
# two spellings of ethanol: 2 unique, 2 duplicates.
if ! diff -w unique.smi "$indir/expected_unique.smi" ||
   ! diff -w duplicate.smi "$indir/expected_duplicate.smi" ; then
  echo "FAIL: unexpected -C results"
  exit 1
fi
echo "PASS"
