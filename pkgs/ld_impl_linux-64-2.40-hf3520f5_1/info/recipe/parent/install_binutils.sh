#!/bin/bash

set -e

cd install

find . -type f -exec bash -c 'mkdir -p /$(dirname {}) && cp {} /{}' ';'

export TARGET="${triplet}"
export OLD_TARGET="${triplet/conda/${ctng_vendor}}"

if [[ "${target_platform}" == win-* ]]; then
  EXEEXT=".exe"
  PREFIX=${PREFIX}/Library
  SYSROOT=${PREFIX}/ucrt64
  OLD_SYSROOT=${PREFIX}/ucrt64
else
  SYSROOT=${PREFIX}/${TARGET}
  OLDSYSROOT=${PREFIX}/${OLD_TARGET}
fi

mkdir -p ${PREFIX}/bin
mkdir -p ${SYSROOT}/bin
mkdir -p ${OLD_SYSROOT}/bin

TOOLS="addr2line ar as c++filt elfedit gprof ld.bfd nm objcopy objdump ranlib readelf size strings strip"

if [[ "${cross_target_platform}" == "linux-"* ]]; then
  TOOLS="${TOOLS} dwp ld.gold"
else
  TOOLS="${TOOLS} dlltool"
fi

# Remove hardlinks and replace them by softlinks
for tool in ${TOOLS}; do
  tool=${tool}${EXEEXT}
  rm -rf ${SYSROOT}/bin/${tool}
  ln -s ${PREFIX}/bin/${TARGET}-${tool} ${SYSROOT}/bin/${tool} || true;
  if [[ "${TARGET}" != "$OLD_TARGET" ]]; then
    ln -s ${PREFIX}/bin/${TARGET}-${tool} ${OLD_SYSROOT}/bin/${tool} || true;
    ln -s ${PREFIX}/bin/${TARGET}-${tool} ${PREFIX}/bin/$OLD_TARGET-${tool} || true;
  fi
  if [[ "$target_platform" == "$cross_target_platform" ]]; then
      mv ${PREFIX}/bin/${tool} ${PREFIX}/bin/${TARGET}-${tool}
  fi
done

rm ${PREFIX}/bin/ld${EXEEXT} || true;
rm ${PREFIX}/bin/${TARGET}-ld${EXEEXT} || true;
rm ${PREFIX}/bin/$OLD_TARGET-ld${EXEEXT} || true;
rm ${OLD_SYSROOT}/bin/ld${EXEEXT} || true;
rm ${SYSROOT}/bin/ld${EXEEXT} || true;

