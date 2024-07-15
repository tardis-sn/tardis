

set -ex



test ! -f ${PREFIX}/lib/libtiff.a
test ! -f ${PREFIX}/lib/libtiffxx.a
test -f ${PREFIX}/lib/libtiff.so
test -f ${PREFIX}/lib/libtiffxx.so
exit 0
