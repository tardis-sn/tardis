

set -ex



cairo-trace --help
test -f $PREFIX/lib/libcairo.a
test -f $PREFIX/lib/libcairo.so
test -f $PREFIX/lib/libcairo-gobject.a
test -f $PREFIX/lib/libcairo-gobject.so
test -f $PREFIX/lib/libcairo-script-interpreter.a
test -f $PREFIX/lib/libcairo-script-interpreter.so
test -f $PREFIX/lib/pkgconfig/cairo-xlib.pc
test -f $PREFIX/include/cairo/cairo.h
grep -q "CAIRO_HAS_FC_FONT 1" $PREFIX/include/cairo/cairo-features.h
export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig
pkg-config --cflags cairo
pkg-config --cflags cairo-fc
pkg-config --cflags cairo-gobject
exit 0
