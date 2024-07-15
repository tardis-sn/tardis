

set -ex



test -f $PREFIX/lib/libedit.so
test ! -f $PREFIX/lib/libedit.a
exit 0
