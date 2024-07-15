#! /bin/bash

set -ex

# Follow the pattern of https://github.com/conda-forge/gdk-pixbuf-feedstock/pull/28/
if [ "${CONDA_BUILD_CROSS_COMPILATION}" = "1" ]; then
    unset _CONDA_PYTHON_SYSCONFIGDATA_NAME
    (
        mkdir -p native-build
        pushd native-build

        export CC=$CC_FOR_BUILD
        export AR=($CC_FOR_BUILD -print-prog-name=ar)
        export NM=($CC_FOR_BUILD -print-prog-name=nm)
        export LDFLAGS=${LDFLAGS//$PREFIX/$BUILD_PREFIX}
        export PKG_CONFIG_PATH=${BUILD_PREFIX}/lib/pkgconfig

        # Unset them as we're ok with builds that are either slow or non-portable
        unset CFLAGS
        unset CPPFLAGS

        meson --libdir=$BUILD_PREFIX/lib --prefix=$BUILD_PREFIX ..
        # This script would generate the functions.txt and dump.xml and save them
        # This is loaded in the native build. We assume that the functions exported
        # by glib are the same for the native and cross builds
        export GI_CROSS_LAUNCHER=$BUILD_PREFIX/libexec/gi-cross-launcher-save.sh
        ninja -j$CPU_COUNT -v
        ninja install
        popd
    )
    export GI_CROSS_LAUNCHER=$BUILD_PREFIX/libexec/gi-cross-launcher-load.sh
fi

# get meson to find pkg-config when cross compiling
export PKG_CONFIG=$BUILD_PREFIX/bin/pkg-config

# need to find gobject-introspection-1.0 as a "native" (build) pkg-config dep
# meson uses PKG_CONFIG_PATH to search when not cross-compiling and
# PKG_CONFIG_PATH_FOR_BUILD when cross-compiling,
# so add the build prefix pkgconfig path to the appropriate variables
export PKG_CONFIG_PATH_FOR_BUILD=$BUILD_PREFIX/lib/pkgconfig
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$BUILD_PREFIX/lib/pkgconfig

meson ${MESON_ARGS:---libdir=$PREFIX/lib} builddir --prefix=$PREFIX
meson configure -Ddocs=false builddir
ninja -v -C builddir
ninja -C builddir install

cd $PREFIX
find . -type f -name "*.la" -exec rm -rf '{}' \; -print
