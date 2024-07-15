

set -ex



test -e $PREFIX/include/blosc.h
test -e $PREFIX/include/blosc-export.h
test ! -e $PREFIX/lib/libblosc.a
test -e $PREFIX/lib/libblosc${SHLIB_EXT}
exit 0
