#!/bin/bash

mkdir build && cd build

cmake ${CMAKE_ARGS} \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_INSTALL_INCLUDEDIR=include \
    ..
make -j${CPU_COUNT}
# remove two tests which need to download data
if [[ "${CONDA_BUILD_CROSS_COMPILATION}" != "1" ]]; then
ctest -E "check_szcomp|sampledata.sh"
fi
make install
