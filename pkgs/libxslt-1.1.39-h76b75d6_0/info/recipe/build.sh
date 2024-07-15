#!/bin/bash

# osx-64 seems to have trouble to get libxml2 info via pkg-config
if [[ ${target_platform} =~ osx.* ]]; then
    LIBXML_CFLAGS="$( pkg-config --cflags libxml-2.0 )"
    LIBXML_LIBS="$( pkg-config --libs libxml-2.0 )"
    export LIBXML_CFLAGS LIBXML_LIBS
fi

./configure --prefix=$PREFIX \
            --with-libxml-prefix=$PREFIX \
            --enable-static=no \
            --without-python

make -j${CPU_COUNT} ${VERBOSE_AT}
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
make check
fi
make install
