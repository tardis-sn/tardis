#!/bin/bash

set -xeo pipefail

if [[ "$target_platform" = osx-* ]] ; then
    # The -dead_strip_dylibs option breaks g-ir-scanner in this package: the
    # scanner links a test executable to find paths to dylibs, but with this
    # option the linker strips them out. The resulting error message is
    # "ERROR: can't resolve libraries to shared libraries: ...".
    export LDFLAGS="$(echo $LDFLAGS |sed -e "s/-Wl,-dead_strip_dylibs//g")"
    export LDFLAGS_LD="$(echo $LDFLAGS_LD |sed -e "s/-dead_strip_dylibs//g")"
fi

# get meson to find pkg-config when cross compiling
export PKG_CONFIG=$BUILD_PREFIX/bin/pkg-config

# need to find gobject-introspection-1.0 as a "native" (build) pkg-config dep
# meson uses PKG_CONFIG_PATH to search when not cross-compiling and
# PKG_CONFIG_PATH_FOR_BUILD when cross-compiling,
# so add the build prefix pkgconfig path to the appropriate variables
export PKG_CONFIG_PATH_FOR_BUILD=$BUILD_PREFIX/lib/pkgconfig
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$BUILD_PREFIX/lib/pkgconfig

export XDG_DATA_DIRS=${XDG_DATA_DIRS}:$PREFIX/share

meson_config_args=(
    -Dintrospection=enabled
    -Dfontconfig=enabled
    -Dfreetype=enabled
    -Dgtk_doc=false
)

# ensure that the post install script is ignored
export DESTDIR="/"

if [[ "$CONDA_BUILD_CROSS_COMPILATION" == "1" ]]; then
  unset _CONDA_PYTHON_SYSCONFIGDATA_NAME
  (
    mkdir -p native-build

    export CC=$CC_FOR_BUILD
    export CXX=$CXX_FOR_BUILD
    export OBJC=$OBJC_FOR_BUILD
    export AR="$($CC_FOR_BUILD -print-prog-name=ar)"
    export NM="$($CC_FOR_BUILD -print-prog-name=nm)"
    export LDFLAGS=${LDFLAGS//$PREFIX/$BUILD_PREFIX}
    export PKG_CONFIG_PATH=${BUILD_PREFIX}/lib/pkgconfig

    # Unset them as we're ok with builds that are either slow or non-portable
    unset CFLAGS
    unset CPPFLAGS
    export host_alias=$build_alias
    export PKG_CONFIG_PATH=$BUILD_PREFIX/lib/pkgconfig

    meson setup native-build \
        "${meson_config_args[@]}" \
        --buildtype=release \
        --prefix=$BUILD_PREFIX \
        -Dlibdir=lib \
        --wrap-mode=nofallback

    # This script would generate the functions.txt and dump.xml and save them
    # This is loaded in the native build. We assume that the functions exported
    # by glib are the same for the native and cross builds
    export GI_CROSS_LAUNCHER=$BUILD_PREFIX/libexec/gi-cross-launcher-save.sh
    ninja -v -C native-build -j ${CPU_COUNT}
    ninja -C native-build install -j ${CPU_COUNT}
  )
  export GI_CROSS_LAUNCHER=$BUILD_PREFIX/libexec/gi-cross-launcher-load.sh
fi

meson setup builddir \
    ${MESON_ARGS} \
    "${meson_config_args[@]}" \
    --buildtype=release \
    --prefix=$PREFIX \
    -Dlibdir=lib \
    --wrap-mode=nofallback
ninja -v -C builddir -j ${CPU_COUNT}
ninja -C builddir install -j ${CPU_COUNT}
