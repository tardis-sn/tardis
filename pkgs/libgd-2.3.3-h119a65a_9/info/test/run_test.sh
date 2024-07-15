

set -ex



test -f $PREFIX/lib/libgd.a
test -f $PREFIX/lib/libgd${SHLIB_EXT}
exit 0
