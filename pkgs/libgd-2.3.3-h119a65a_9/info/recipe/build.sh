#!/usr/bin/env bash

set -x

export TMPDIR=$SRC_DIR/temp_files
mkdir -p $TMPDIR

if [[ ${target_platform} != *-64 ]]; then
  # https://github.com/libgd/libgd/issues/278
  export CFLAGS="$CFLAGS -ffp-contract=off"
elif [[ ${target_platform} == linux-32 ]]; then
  # https://github.com/libgd/libgd/blob/gd-2.2.5/docs/README.TESTING#L65-L70
  export CFLAGS="$CFLAGS -msse -mfpmath=sse"
fi

autoreconf -vfi
./configure --prefix=$PREFIX \
            --without-xpm \
            --without-x \
            --disable-werror \
|| { cat config.log; exit 1; }

make -j${CPU_COUNT} ${VERBOSE_AT}
make install

if [[ "${CONDA_BUILD_CROSS_COMPILATION}" != "1" ]]; then
  if [[ $target_platform == osx-64 ]]; then
    make check || failed=1
  grep -rl "DYLD_LIBRARY_PATH=" tests | xargs sed -i.bak "s~DYLD_LIBRARY_PATH=.*~DYLD_LIBRARY_PATH=$PREFIX/lib~g"
  fi

  #^ see: https://github.com/libgd/libgd/issues/302
  export FREETYPE_PROPERTIES=truetype:interpreter-version=35
  make check || { cat tests/test-suite.log; exit 1; }
fi

# We can remove this when we start using the new conda-build.
find $PREFIX -name '*.la' -delete
