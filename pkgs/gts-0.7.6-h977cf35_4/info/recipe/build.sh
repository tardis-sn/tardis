#!/bin/sh
# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/gnuconfig/config.* .

./configure --prefix=$PREFIX --with-pic

make -j$CPU_COUNT
make install
