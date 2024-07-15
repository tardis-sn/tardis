

set -ex



test -e $PREFIX/include/zlib-ng.h
test -e $PREFIX/include/zconf-ng.h
test -e $PREFIX/lib/libz-ng${SHLIB_EXT}
test -e $PREFIX/lib/pkgconfig/zlib-ng.pc
exit 0
