#!/bin/bash

set -ex

# necessary to ensure the gobject-introspection-1.0 pkg-config file gets found
# meson needs this to determine where the g-ir-scanner script is located
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$BUILD_PREFIX/lib/pkgconfig
export PKG_CONFIG=$BUILD_PREFIX/bin/pkg-config
# Make sure .gir files in $PREFIX are found
export XDG_DATA_DIRS=${XDG_DATA_DIRS}:$PREFIX/share:$BUILD_PREFIX/share

meson_config_args=(
    -Dbenchmark=disabled
    -Dcairo=enabled
    -Dcoretext=auto
    -Ddocs=disabled
    -Dfreetype=enabled
    -Dgdi=auto
    -Dglib=enabled
    -Dgobject=enabled
    -Dgraphite2=enabled
    -Dicu=enabled
    -Dintrospection=enabled
    -Dtests=disabled
)

if [ -n "$OSX_ARCH" ] ; then
    # The -dead_strip_dylibs option breaks g-ir-scanner here
    export LDFLAGS="$(echo $LDFLAGS |sed -e "s/-Wl,-dead_strip_dylibs//g")"
    export LDFLAGS_LD="$(echo $LDFLAGS_LD |sed -e "s/-dead_strip_dylibs//g")"
fi

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

    meson setup native-build \
        "${meson_config_args[@]}" \
        --buildtype=release \
        --prefix=$BUILD_PREFIX \
        -Dlibdir=lib \
        --wrap-mode=nofallback

    # This script would generate the functions.txt and dump.xml and save them
    # This is loaded in the native build. We assume that the functions exported
    # by the package are the same for the native and cross builds
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
    --default-library=both \
    --prefix=$PREFIX \
    -Dlibdir=lib \
    --wrap-mode=nofallback
ninja -v -C builddir -j ${CPU_COUNT}
ninja -C builddir install -j ${CPU_COUNT}


pushd $PREFIX
rm -rf share/gtk-doc
popd
