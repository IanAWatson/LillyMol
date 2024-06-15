#!/bin/bash
# Install LillyMol executables, shell wrappers, data and shared libraries
# to LILLYMOL_HOME, which must be defined.

if [[ ! -v LILLYMOL_HOME ]] ; then
  echo "Shell variable LILLYMOL_HOME must be defined" >&2
  exit 1
fi

bindir="${LILLYMOL_HOME}/bin/$(uname)"
if [[ ! -d "${bindir}" ]] ; then
  mkdir -p ${bindir}
fi

for src in bin/$(uname)/* ; do
  name=$(basename $src)
  dest="${bindir}/${name}"
  if [[ ! -s ${dest} || ${src} -nt ${dest} ]] ; then
    echo ${src}
    cp -f $src $dest
  fi
done

# Special processing for outside Lilly and svmfp models.
if [[ $(hostname -d) =~ 'lilly.com' ]] ; then
else
  cp ${LILLYMOL_HOME}/bin/$(uname)/gfp_to_svm_lite_v3 ${LILLYMOL_HOME}/bin/$(uname)/gfp_to_svm_lite
fi

echo 'Scripts' >&2
for src in contrib/bin/*.{sh,rb,py,pl} ; do
  name=$(basename $src)
  dest="${LILLYMOL_HOME}/bin/${name}"
  if [[ ! -s ${dest} || ${src} -nt ${dest} ]] ; then
    cp -f $src $dest
    echo $name
  fi
done

echo 'Data' >&2
cp -r data ${LILLYMOL_HOME}

echo 'Lib' >@2
cp -f bin/lib ${LILLYMOL_HOME}/bin

libdir="${LILLYMOL_HOME}/lib"
if [[ ! -d "${libdir}" ]] ; then
  mkdir -p "${libdir}"
fi

echo "Shared libraries"
for src in lib/*.so ; do
  name=$(basename $src)
  dest="${LILLYMOL_HOME}/lib/${name}"
  if [[ ! -s ${dest} || ${src} -nt ${dest} ]] ; then
    cp -f $src $dest
    echo $name
  fi
done

echo "Compiled Protos" >&2

if [[ ! -d "${LILLYMOL_HOME}/proto" ]] ; then
  mkdir -p "${LILLYMOL_HOME}/proto"
fi

for src in $(find src -name '*_pb2.py') ; do
  name=$(basename $src)
  dest="${LILLYMOL_HOME}/proto/${name}"
  if [[ ! -s ${dest} || ${src} -nt ${dest} ]] ; then
    cp -f $src $dest
    echo $name
  fi
done
