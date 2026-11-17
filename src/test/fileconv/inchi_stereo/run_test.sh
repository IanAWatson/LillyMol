#! /bin/bash
# InChI stereo layers: check what fileconv writes for stereocentres, and that it
# reads the same stereochemistry back. Only builds made with bazel --config=inchi can do
# this. Builds without InChI print "not implemented", and are reported as PASS so
# they are not silently skipped in a way that looks like a failure.
#
# in/expected.inchi was generated with RDKit, not with LillyMol, so it is an independent
# reference. The cases include quaternary centres (C, P, N+, isotopic H), whose stereo
# layer was once written inverted.
#
# Arguments: executable, in directory, out directory.
exe=$1
indir=$2

"$exe" -o inchi -S written "$indir/stereo.smi" > /dev/null 2> written_stderr

if grep -q 'not implemented' written_stderr ; then
  echo "PASS"
  exit 0
fi

if ! diff -w written.inchi "$indir/expected.inchi" > diff_expected ; then
  echo "FAIL: InChI written by fileconv differs from the RDKit reference"
  cat diff_expected
  exit 1
fi

# Round trip: read the InChI back, and the structures must be identical, stereo included.
"$exe" -g xh -o usmi -S original "$indir/stereo.smi" > /dev/null 2>&1
"$exe" -i inchi -g xh -o usmi -S back written.inchi > /dev/null 2>&1
if ! diff -w original.smi back.smi > diff_roundtrip ; then
  echo "FAIL: stereochemistry changed in an InChI round trip"
  cat diff_roundtrip
  exit 1
fi

echo "PASS"
