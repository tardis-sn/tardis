

set -ex



brotli --help
test -f $PREFIX/lib/libbrotlienc$SHLIB_EXT
test -f $PREFIX/lib/libbrotlidec$SHLIB_EXT
test -f $PREFIX/lib/libbrotlicommon$SHLIB_EXT
test -f $PREFIX/include/brotli/encode.h
test ! -f $PREFIX/lib/libbrotlienc-static.a
test ! -f $PREFIX/lib/libbrotlidec-static.a
test ! -f $PREFIX/lib/libbrotlicommon-static.a
exit 0
