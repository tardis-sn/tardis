#!/bin/bash
# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/libtool/build-aux/config.* ./config

# Pass explicit paths to the prefix for each dependency.
./configure --prefix="${PREFIX}" \
            --disable-static \
            --host=$HOST \
            --build=$BUILD \
            --with-zlib-include-dir="${PREFIX}/include" \
            --with-zlib-lib-dir="${PREFIX}/lib" \
            --with-jpeg-include-dir="${PREFIX}/include" \
            --with-jpeg-lib-dir="${PREFIX}/lib" \
            --with-lzma-include-dir="${PREFIX}/include" \
            --with-lzma-lib-dir="${PREFIX}/lib" \
            --with-zstd-include-dir="${PREFIX}/include" \
            --with-zstd-lib-dir="${PREFIX}/lib"

make -j${CPU_COUNT} ${VERBOSE_AT}
if [[ "$CONDA_BUILD_CROSS_COMPILATION" != 1 ]]; then
  make check
fi
make install

rm -rf "${TIFF_BIN}" "${TIFF_SHARE}" "${TIFF_DOC}"

# For some reason --docdir is not respected above.
rm -rf "${PREFIX}/share"

# We can remove this when we start using the new conda-build.
find $PREFIX -name '*.la' -delete
