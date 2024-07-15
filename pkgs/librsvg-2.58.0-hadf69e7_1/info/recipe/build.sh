#! /bin/bash
set -ex

# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/gnuconfig/config.* .

# $BUILD_PREFIX needed here so gi-docgen can find .gir files:
export XDG_DATA_DIRS=${XDG_DATA_DIRS}:$PREFIX/share:$BUILD_PREFIX/share

# https://github.com/rust-lang/cargo/issues/10583#issuecomment-1129997984
export CARGO_NET_GIT_FETCH_WITH_CLI=true

configure_args=(
    --disable-Bsymbolic
    --disable-static
    --enable-pixbuf-loader=yes
    --enable-introspection=yes
)

if [[ $target_platform == osx-* ]] ; then
  # Workaround for https://gitlab.gnome.org/GNOME/librsvg/-/issues/545 ; should be removable soon.
  export LDFLAGS="$LDFLAGS -lobjc"
fi


if [[ "$CONDA_BUILD_CROSS_COMPILATION" == 1 ]]; then
  unset _CONDA_PYTHON_SYSCONFIGDATA_NAME
  (
    unset CARGO_BUILD_TARGET
    mkdir -p native-build
    pushd native-build

    export CC=$CC_FOR_BUILD
    export AR="$($CC_FOR_BUILD -print-prog-name=ar)"
    export NM="$($CC_FOR_BUILD -print-prog-name=nm)"
    export LDFLAGS=${LDFLAGS//$PREFIX/$BUILD_PREFIX}
    export PKG_CONFIG_PATH=${BUILD_PREFIX}/lib/pkgconfig

    # Unset them as we're ok with builds that are either slow or non-portable
    unset CFLAGS
    unset CPPFLAGS
    export host_alias=$build_alias
    export PKG_CONFIG_PATH=$BUILD_PREFIX/lib/pkgconfig

    ../configure --prefix=$BUILD_PREFIX "${configure_args[@]}"

    # This script would generate the functions.txt and dump.xml and save them
    # This is loaded in the native build. We assume that the functions exported
    # by glib are the same for the native and cross builds
    export GI_CROSS_LAUNCHER=$BUILD_PREFIX/libexec/gi-cross-launcher-save.sh
    make -j${CPU_COUNT}
    make install
    popd
  )
  export GI_CROSS_LAUNCHER=$BUILD_PREFIX/libexec/gi-cross-launcher-load.sh
fi

export PKG_CONFIG_PATH_FOR_BUILD=$BUILD_PREFIX/lib/pkgconfig
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$BUILD_PREFIX/lib/pkgconfig

export RUST_TARGET=$CARGO_BUILD_TARGET
unset CARGO_BUILD_TARGET

./configure --prefix=$PREFIX "${configure_args[@]}" || { cat config.log ; exit 1 ; }
make -j$CPU_COUNT
make install

rm -rf $PREFIX/share/gtk-doc
