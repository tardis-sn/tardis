

set -ex



test -e $PREFIX/include/snappy.h
test -e $PREFIX/include/snappy-stubs-public.h
test -e $PREFIX/lib/libsnappy$SHLIB_EXT
test -f $PREFIX/lib/libsnappy.so.${PKG_VERSION}
exit 0
