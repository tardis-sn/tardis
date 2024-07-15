#!/bin/bash

set -xeo pipefail

# For the Linux builds, get our prefix into the shared library
sed -i.orig s:'@PREFIX@':"$PREFIX":g src/fccfg.c

# osx-arm64 cross build sets up $MESON_ARGS with these settings plus the cross file:
if [[ $MESON_ARGS == "" ]] ; then
    MESON_ARGS="--buildtype=release --prefix=$PREFIX -Dlibdir=lib"
fi

if [[ "$target_platform" == "osx-arm64" && "$CONDA_BUILD_CROSS_COMPILATION" == "1" ]]; then
    export PKG_CONFIG=$BUILD_PREFIX/bin/pkg-config

    # When cross-compiling the environment is set up as if Python will be used
    # as a host/run tool, but for us it is genuinely just a build tool. This
    # override fixes Meson's detection of Python given that:
    export _CONDA_PYTHON_SYSCONFIGDATA_NAME=_sysconfigdata_x86_64_apple_darwin13_4_0
fi

meson_setup_args=(
    $MESON_ARGS
    --default-library=both
    --wrap-mode=nofallback
)

meson setup builddir "${meson_setup_args[@]}" || { cat builddir/meson-logs/meson-log.txt ; exit 1 ; }

ninja -v -C builddir -j ${CPU_COUNT}

if [[ "$CONDA_BUILD_CROSS_COMPILATION" != "1" ]]; then
    ninja -C builddir test -j ${CPU_COUNT}
fi

ninja -C builddir install -j ${CPU_COUNT}

# Clear out the local cache but make sure the directory is packaged.
rm -Rf "$PREFIX"/var/cache/fontconfig
mkdir -p "$PREFIX"/var/cache/fontconfig
touch "$PREFIX"/var/cache/fontconfig/.leave
