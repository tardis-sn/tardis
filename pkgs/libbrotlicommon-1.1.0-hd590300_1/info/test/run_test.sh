

set -ex



test -f $PREFIX/lib/libbrotlicommon$SHLIB_EXT
test ! -f $PREFIX/lib/libbrotlicommon-static.a
exit 0
