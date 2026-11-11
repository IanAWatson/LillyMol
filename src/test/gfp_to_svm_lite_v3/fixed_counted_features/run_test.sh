#!/bin/bash
# Fixed size counted fingerprints are converted to svm_lite features, just like
# sparse counted ones. Feature numbers depend on hash ordering, so write each row
# as molecule, activity, and the sorted counts, which do not.
# $1 is the executable, $2 the 'in' directory. Runs in the temporary directory.

exe=$1
indir=$2

canonical() {
  while IFS= read -r line; do
    body=${line%%#*}
    name=${line#*# }
    set -- ${body}
    activity=$1
    shift
    counts=$(for token in "$@"; do echo "${token#*:}"; done | sort -n | tr '\n' ' ')
    echo "${name} ${activity} ${counts}"
  done < "$1"
}

echo "fixed counted"
"${exe}" -A "${indir}/act.txt" -O svml -C fc -S fc "${indir}/in.fc.gfp" 2> /dev/null || echo "FAILED"
canonical fc.svml

echo "sparse counted"
"${exe}" -A "${indir}/act.txt" -O svml -C nc -S nc "${indir}/in.nc.gfp" 2> /dev/null || echo "FAILED"
canonical nc.svml

if cmp -s <(canonical fc.svml) <(canonical nc.svml); then
  echo "fixed counted and sparse counted agree"
else
  echo "fixed counted and sparse counted differ"
fi

# Features are known by fingerprint tag, so a cross reference made from the sparse
# form cannot be applied to the fixed counted form. This must fail, not write
# records with no features.
if "${exe}" -r -O svml -X nc_xref.dat -S mismatch "${indir}/in.fc.gfp" > /dev/null 2>&1; then
  echo "mismatched cross reference: unexpected success"
else
  echo "mismatched cross reference: failed as expected"
fi

# A cross reference from the fixed counted form can be applied to it.
"${exe}" -r -O svml -X fc_xref.dat -S applied "${indir}/in.fc.gfp" 2> /dev/null || echo "FAILED"
echo "applied rows $(wc -l < applied.svml)"

# Every error must give a non-zero exit status.
expect_failure() {
  description=$1
  shift
  if "$@" > /dev/null 2>&1; then
    echo "${description}: unexpected success"
  else
    echo "${description}: failed as expected"
  fi
}

expect_failure "unwritable output" "${exe}" -r -O svml -X fc_xref.dat -S /no/such/dir/out "${indir}/in.fc.gfp"
expect_failure "unwritable proto" "${exe}" -r -O svml -C /no/such/dir/proto -S x "${indir}/in.fc.gfp"
expect_failure "gfp output not implemented" "${exe}" -r -O gfp -C g -S g "${indir}/in.fc.gfp"
expect_failure "support removes everything" "${exe}" -r -O svml -p 100 -C h -S h "${indir}/in.fc.gfp"
expect_failure "empty input" "${exe}" -r -O svml -C e -S e /dev/null

# The input is read twice when creating a cross reference.
cat "${indir}/in.fc.gfp" > stdin_copy.gfp
expect_failure "cross reference from stdin" bash -c "cat stdin_copy.gfp | '${exe}' -r -O svml -C s -S s -"

# A record with no identifier must not be silently dropped.
sed '0,/^PCN</s//XXX</' "${indir}/in.fc.gfp" > noid.gfp
expect_failure "record with no identifier" "${exe}" -r -O svml -C n -S n noid.gfp

# Applying a cross reference made from two fingerprints to an input that has
# only one of them must fail, not write records missing those features.
awk 'NR==FNR { if ($0 ~ /^FCFML</) fc[++n] = $0; next }
     { if ($0 == "|") { m++; print fc[m] } print }' "${indir}/in.fc.gfp" "${indir}/in.nc.gfp" > both.gfp
"${exe}" -A "${indir}/act.txt" -O svml -C both -S both both.gfp 2> /dev/null || echo "FAILED"
expect_failure "cross reference has an extra fingerprint" "${exe}" -r -O svml -X both_xref.dat -S extra "${indir}/in.fc.gfp"

# Applying a cross reference to input read from stdin works.
cat "${indir}/in.fc.gfp" | "${exe}" -r -O svml -X fc_xref.dat -S from_stdin - 2> /dev/null || echo "FAILED"
echo "applied from stdin rows $(wc -l < from_stdin.svml)"

# A feature is kept if it is found in at least -p molecules. The sulphur and the
# halogens are each only in one molecule, so -p 2 removes them.
nfeatures() {
  while IFS= read -r line; do
    body=${line%%#*}
    name=${line#*# }
    set -- ${body}
    shift
    echo "${name} $#"
  done < "$1"
}
"${exe}" -A "${indir}/act.txt" -O svml -p 1 -C p1 -S p1 "${indir}/in.fc.gfp" 2> /dev/null || echo "FAILED"
"${exe}" -A "${indir}/act.txt" -O svml -p 2 -C p2 -S p2 "${indir}/in.fc.gfp" 2> /dev/null || echo "FAILED"
echo "features per molecule, support 1 then 2"
paste -d ' ' <(nfeatures p1.svml) <(nfeatures p2.svml | cut -d ' ' -f 2)

# Integer molecular properties (MPR). Each non zero property is a feature whose
# value is the property value. They are not subject to the support requirement,
# so with a high support level only the properties remain.
echo "with molecular properties"
"${exe}" -A "${indir}/act.txt" -O svml -C mpr -S mpr "${indir}/in.mpr.gfp" 2> /dev/null || echo "FAILED"
canonical mpr.svml
"${exe}" -A "${indir}/act.txt" -O svml -p 100 -C mprp -S mprp "${indir}/in.mpr.gfp" 2> /dev/null || echo "FAILED"
echo "with molecular properties only, support 100"
canonical mprp.svml

# The cross reference and the input must agree about the properties.
expect_failure "cross reference with properties, input without" "${exe}" -r -O svml -X mpr_xref.dat -S m1 "${indir}/in.fc.gfp"
expect_failure "cross reference without properties, input with" "${exe}" -r -O svml -X fc_xref.dat -S m2 "${indir}/in.mpr.gfp"
"${exe}" -r -O svml -X mpr_xref.dat -S m3 "${indir}/in.mpr.gfp" 2> /dev/null || echo "FAILED"
echo "applied with properties rows $(wc -l < m3.svml)"
