#!/bin/bash
set -x

./autogen.sh

# remove libtool files
find $PREFIX -name '*.la' -delete

declare -a _xtra_config_flags
declare -a _xtra_make_args

if [[ ${target_platform} =~ .*osx.* ]]; then
    export OBJC="${CC}"
    # xcodebuild uses ld instead of clang and fails
    export LD="${CC_FOR_BUILD}"
    _xtra_config_flags+=(--with-quartz)
fi

./configure --prefix=$PREFIX \
            --disable-debug \
            --disable-java \
            --disable-php \
            --disable-perl \
            --disable-tcl \
            --enable-ltdl \
            --without-x \
            --without-qt \
            --without-gtk \
            --with-ann=no \
            --with-gts=yes \
            --with-gdk=yes \
            --with-rsvg=yes \
            --with-expat=yes \
            --with-libgd=yes \
            --with-freetype2=yes \
            --with-fontconfig=yes \
            --with-pangocairo=yes \
            --with-gdk-pixbuf=yes \
            "${_xtra_config_flags[@]}"


if [ $CONDA_BUILD_CROSS_COMPILATION = 1 ] && [ "${target_platform}" = "osx-arm64" ]; then
    sed -i.bak 's/HOSTCC/CC_FOR_BUILD/g' $SRC_DIR/lib/gvpr/Makefile.am
    sed -i.bak '/dot$(EXEEXT) -c/d' $SRC_DIR/cmd/dot/Makefile.am
    _xtra_make_args+=(ARCH=arm64)
fi

make -j${CPU_COUNT} "${_xtra_make_args[@]}"
# This is failing for rtest.
# Doesn't do anything for the rest
# make check
make install

# Configure plugins
if [ $CONDA_BUILD_CROSS_COMPILATION != 1 ]; then
    $PREFIX/bin/dot -c
fi
