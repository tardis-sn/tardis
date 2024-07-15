#!/bin/bash

rm -rf internal-complibs

mkdir build
cd build

cmake ${CMAKE_ARGS} -G "Unix Makefiles" \
      -DCMAKE_BUILD_TYPE="Release" \
      -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
      -DCMAKE_POSITION_INDEPENDENT_CODE=1 \
      -DBUILD_STATIC=0 \
      -DBUILD_SHARED=1 \
      -DBUILD_TESTS=1 \
      -DBUILD_BENCHMARKS=0 \
      -DDEACTIVATE_SNAPPY=0 \
      -DPREFER_EXTERNAL_LZ4=ON \
      -DPREFER_EXTERNAL_ZLIB=ON \
      -DPREFER_EXTERNAL_ZSTD=ON \
      "${SRC_DIR}"

cmake --build .
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
ctest
fi
cmake --build . --target install
