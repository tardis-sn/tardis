#!/bin/bash
set -ex

# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/gnuconfig/config.* .

# The libwebp build script doesn't pick all the other libraries up on its own
# (even though it should by using PREFIX), so pass all the necessary parameters
# for finding other imaging libraries to the configure script.
./configure \
    --disable-dependency-tracking \
    --disable-gl \
    --disable-static \
    --enable-libwebpdecoder \
    --enable-libwebpdemux \
    --enable-libwebpmux \
    --prefix=${PREFIX} \
    --with-gifincludedir=${PREFIX}/include \
    --with-giflibdir=${PREFIX}/lib \
    --with-jpegincludedir=${PREFIX}/include \
    --with-jpeglibdir=${PREFIX}/lib \
    --with-tiffincludedir=${PREFIX}/include \
    --with-tifflibdir=${PREFIX}/lib \

make

if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
make check
fi

make install
