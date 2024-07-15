

set -ex



test -e $PREFIX/include/blosc2.h
test -e $PREFIX/include/blosc2/blosc2-export.h
test -e $PREFIX/lib/libblosc2${SHLIB_EXT}
test -e $PREFIX/lib/pkgconfig/blosc2.pc
exit 0
