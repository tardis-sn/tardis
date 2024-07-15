#!/bin/bash

source ${RECIPE_DIR}/setup_compiler.sh
set -e -x

export CHOST="${triplet}"

# libtool wants to use ranlib that is here, macOS install doesn't grok -t etc
# .. do we need this scoped over the whole file though?
#export PATH=${SRC_DIR}/gcc_built/bin:${SRC_DIR}/.build/${CHOST}/buildtools/bin:${SRC_DIR}/.build/tools/bin:${PATH}

pushd ${SRC_DIR}/build

make -C ${CHOST}/libgcc prefix=${PREFIX} install

# ${PREFIX}/lib/libgcc_s.so* goes into libgcc-ng output, but
# avoid that the equivalents in ${PREFIX}/${CHOST}/lib end up
# in gcc_impl_{{ cross_target_platform }}, c.f. install-gcc.sh
mkdir -p ${PREFIX}/${CHOST}/lib
mv ${PREFIX}/lib/libgcc_s.so* ${PREFIX}/${CHOST}/lib
# This is in gcc_impl as it is gcc specific and clang has the same header
rm -rf ${PREFIX}/lib/gcc/${CHOST}/${gcc_version}/include/unwind.h

popd

