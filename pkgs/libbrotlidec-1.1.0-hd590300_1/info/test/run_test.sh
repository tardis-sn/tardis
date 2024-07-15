

set -ex



test -f $PREFIX/lib/libbrotlidec$SHLIB_EXT
test ! -f $PREFIX/lib/libbrotlidec-static.a
exit 0
