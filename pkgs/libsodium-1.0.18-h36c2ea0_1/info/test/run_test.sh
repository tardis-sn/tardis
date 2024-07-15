

set -ex



test -f ${PREFIX}/include/sodium.h
test -f ${PREFIX}/lib/libsodium.a
test -f ${PREFIX}/lib/libsodium.so
exit 0
