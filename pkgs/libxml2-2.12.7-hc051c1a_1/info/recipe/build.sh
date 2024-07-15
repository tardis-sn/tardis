#!/bin/bash
# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/libtool/build-aux/config.* .

./autogen.sh

if [[ `uname` == "Linux" ]]; then
  # workaround weird configure behaviour where it decides
  # it doesn't need libiconv
  export LDFLAGS="-L${PREFIX}/lib -liconv"
fi

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$PREFIX/lib/pkgconfig

./configure --prefix="${PREFIX}" \
            --build=${BUILD} \
            --host=${HOST} \
            --with-iconv="${PREFIX}" \
            --with-zlib="${PREFIX}" \
            --with-icu \
            --with-lzma="${PREFIX}" \
            --with-ftp \
            --with-legacy \
            --with-python=no \
            --with-tls \
            --enable-static=no || cat config.log
make -j${CPU_COUNT} ${VERBOSE_AT}

if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
  make check ${VERBOSE_AT}
fi

make install

# Remove large documentation files that can take up to 6.6/9.2MB of the install
# size.
# https://github.com/conda-forge/libxml2-feedstock/issues/57
rm -rf ${PREFIX}/share/doc
rm -rf ${PREFIX}/share/gtk-doc
rm -rf ${PREFIX}/share/man

for f in "activate" "deactivate"; do
    mkdir -p "${PREFIX}/etc/conda/${f}.d"
    cp "${RECIPE_DIR}/${f}.sh" "${PREFIX}/etc/conda/${f}.d/${PKG_NAME}_${f}.sh"
done
