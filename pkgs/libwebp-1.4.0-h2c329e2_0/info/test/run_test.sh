

set -ex



test -f $PREFIX/bin/cwebp
test -f $PREFIX/bin/dwebp
test -f $PREFIX/bin/gif2webp
test -f $PREFIX/bin/img2webp
test -f $PREFIX/bin/webpinfo
test -f $PREFIX/bin/webpmux
test ! -f $PREFIX/lib/libwebdecoder.a
test ! -f $PREFIX/lib/libwebp.a
test ! -f $PREFIX/lib/libwebpdemux.a
test ! -f $PREFIX/lib/libwebpmux.a
dwebp examples/test.webp
exit 0
