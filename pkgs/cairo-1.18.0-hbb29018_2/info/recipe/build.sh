#!/bin/bash

set -ex

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$BUILD_PREFIX/lib/pkgconfig
export PKG_CONFIG=$BUILD_PREFIX/bin/pkg-config

meson_config_args=(
    -Dfontconfig=enabled
    -Dfreetype=enabled
    -Dglib=enabled
)

if test $(uname) == Darwin ; then
    meson_config_args+=(
        -Dxlib=disabled
        -Dxlib-xcb=disabled
        -Dxcb=disabled
    )
elif test $(uname) == Linux ; then
    meson_config_args+=(-Dxlib-xcb=enabled)
fi

if [[ "$CONDA_BUILD_CROSS_COMPILATION" == 1 ]]; then
    # See: https://gitlab.freedesktop.org/cairo/cairo/-/merge_requests/134
    cat <<EOF >cross_file.txt
[properties]
ipc_rmid_deferred_release = true
EOF
    meson_config_args+=(--cross-file cross_file.txt)
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
