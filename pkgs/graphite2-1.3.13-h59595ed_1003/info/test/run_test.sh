

set -ex



test -f $PREFIX/lib/pkgconfig/graphite2.pc
test -f $PREFIX/share/graphite2/graphite2.cmake
test -f $PREFIX/include/graphite2/Font.h
test -f $PREFIX/lib/libgraphite2${SHLIB_EXT}
conda inspect linkages -p $PREFIX $PKG_NAME
exit 0
