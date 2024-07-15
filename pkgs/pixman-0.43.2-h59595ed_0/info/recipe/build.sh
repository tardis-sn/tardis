#!/bin/bash

set -ex

meson_config_args=()

if [[ $(uname) == Darwin ]]; then
  meson_config_args+=(-Dopenmp=disabled)
fi

if [ "${target_platform}" == linux-ppc64le ]; then
  meson_config_args+=(-Dvmx=disabled)
fi

if [ "${target_platform}" == osx-arm64 ]; then
  meson_config_args+=(-Da64-neon=disabled)
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
