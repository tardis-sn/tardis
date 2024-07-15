

set -ex



test -f $PREFIX/lib/libbrotlienc$SHLIB_EXT
test ! -f $PREFIX/lib/libbrotlienc-static.a
exit 0
