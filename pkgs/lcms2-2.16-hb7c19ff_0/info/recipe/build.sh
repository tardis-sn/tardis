#!/bin/bash
# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/libtool/build-aux/config.* .

./configure \
    --prefix=$PREFIX \
    --disable-static \
    --with-tiff=$PREFIX \
    --with-jpeg=$PREFIX

make -j${CPU_COUNT}
# 202210 - hmaarrfk
# Without installing before make check, it seemed that
# make check would fail on OSX-64bit
make install
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
make check
fi
