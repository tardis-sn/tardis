#!/bin/bash

declare -a CMAKE_PLATFORM_FLAGS
if [[ ${target_platform} == osx-64 ]]; then
  CMAKE_PLATFORM_FLAGS+=(-DCMAKE_OSX_SYSROOT="${CONDA_BUILD_SYSROOT}")
fi

mkdir build || true
pushd build

  cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DTIFF_LIBRARY=$PREFIX/lib/libtiff${SHLIB_EXT} \
        -DTIFF_INCLUDE_DIR=$PREFIX/include \
        -DPNG_LIBRARY_RELEASE=$PREFIX/lib/libpng${SHLIB_EXT} \
        -DPNG_PNG_INCLUDE_DIR=$PREFIX/include \
        -DZLIB_LIBRARY=$PREFIX/lib/libz${SHLIB_EXT} \
        -DZLIB_INCLUDE_DIR=$PREFIX/include \
        -DBUILD_JPWL=OFF \
        "${CMAKE_PLATFORM_FLAGS[@]}" \
        $SRC_DIR

  make -j${CPU_COUNT} ${VERBOSE_CM}
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
  ctest
fi
  make install -j${CPU_COUNT}

popd