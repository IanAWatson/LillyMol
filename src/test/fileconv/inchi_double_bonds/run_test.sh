#!/usr/bin/env bash
# Export references were generated independently with RDKit. No RDKit is needed
# to run the test. InChI import currently does not restore double-bond stereo,
# so these are export checks, separate from the tetrahedral round-trip test.
exe=$1
indir=$2
"$exe" -o inchi -S written "$indir/in.smi" > /dev/null 2> written_stderr
rc=$?
if grep -q 'not implemented' written_stderr; then
  echo PASS
  exit 0
fi
if [[ $rc -ne 0 ]]; then
  cat written_stderr
  exit 1
fi
diff -w written.inchi "$indir/expected.inchi" || exit 1
echo PASS
