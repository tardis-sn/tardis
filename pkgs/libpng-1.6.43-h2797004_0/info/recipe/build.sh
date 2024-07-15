#!/bin/bash
# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/libtool/build-aux/config.* .

export CFLAGS="$CFLAGS -I$PREFIX/include -L$PREFIX/lib"
export CPPFLAGS="$CPPFLAGS -I$PREFIX/include"

./configure --prefix=$PREFIX \
            --with-zlib-prefix=$PREFIX

make -j${CPU_COUNT} ${VERBOSE_AT}
if [[ "${CONDA_BUILD_CROSS_COMPILATION}" != "1" ]]; then
make check
fi
make install

# We can remove this when we start using the new conda-build.
find $PREFIX -name '*.la' -delete
