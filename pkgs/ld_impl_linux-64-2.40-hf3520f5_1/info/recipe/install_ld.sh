#!/bin/bash

set -e

TARGET="${triplet}"
OLD_TARGET="${triplet/conda/${ctng_vendor}}"

if [[ "$target_platform" == win-* ]]; then
  EXEEXT=".exe"
  PREFIX=$PREFIX/Library
  SYSROOT=$PREFIX/ucrt64
  OLD_SYSROOT=$PREFIX/ucrt64
else
  SYSROOT=$PREFIX/${TARGET}
  OLD_SYSROOT=$PREFIX/${OLD_TARGET}
fi

mkdir -p $PREFIX/bin
mkdir -p $OLD_SYSROOT/bin
mkdir -p $SYSROOT/bin

if [[ "$target_platform" == "$cross_target_platform" ]]; then
  cp $PWD/install/$PREFIX/bin/ld${EXEEXT} $PREFIX/bin/$TARGET-ld${EXEEXT}
else
  cp $PWD/install/$PREFIX/bin/$TARGET-ld${EXEEXT} $PREFIX/bin/$TARGET-ld${EXEEXT}
fi

if [[ "$TARGET" != "$OLD_TARGET" ]]; then
  ln -s $PREFIX/bin/$TARGET-ld${EXEEXT} $PREFIX/bin/$OLD_TARGET-ld${EXEEXT}
  ln -s $PREFIX/bin/$TARGET-ld${EXEEXT} $OLD_SYSROOT/bin/ld${EXEEXT}
fi
ln -s $PREFIX/bin/$TARGET-ld${EXEEXT} $SYSROOT/bin/ld${EXEEXT}
