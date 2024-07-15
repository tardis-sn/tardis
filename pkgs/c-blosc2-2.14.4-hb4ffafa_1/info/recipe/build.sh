#!/bin/bash

rm -rf internal-complibs

mkdir build
cd build

if [[ "${target_platform}" != "${build_platform}" ]]; then
    BUILD_TESTS=0
else
    BUILD_TESTS=1
fi

cmake -G "Unix Makefiles" \
      ${CMAKE_ARGS} \
      -DCMAKE_BUILD_TYPE="Release" \
      -DCMAKE_POSITION_INDEPENDENT_CODE=1 \
      -DBUILD_STATIC=0 \
      -DBUILD_SHARED=1 \
      -DBUILD_TESTS=${BUILD_TESTS} \
      -DBUILD_EXAMPLES=0 \
      -DBUILD_BENCHMARKS=0 \
      -DPREFER_EXTERNAL_LZ4=1 \
      -DPREFER_EXTERNAL_ZSTD=1 \
      -DPREFER_EXTERNAL_ZLIB=1 \
      "${SRC_DIR}"

cmake --build .
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
    ctest
fi
cmake --build . --target install
