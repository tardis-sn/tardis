

set -ex



test ! -f ${PREFIX}/lib/libdeflate.a
test -f ${PREFIX}/lib/libdeflate${SHLIB_EXT}
test -f ${PREFIX}/include/libdeflate.h
libdeflate-gzip -h
libdeflate-gunzip -h
exit 0
