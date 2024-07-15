

set -ex



gdk-pixbuf-csource --version
gdk-pixbuf-query-loaders
gdk-pixbuf-pixdata --help
test -f $PREFIX/lib/libgdk_pixbuf-2.0${SHLIB_EXT}
test -f `pkg-config --variable=libdir --dont-define-prefix gdk-pixbuf-2.0`/libgdk_pixbuf-2.0${SHLIB_EXT}
test -f $PREFIX/lib/girepository-1.0/GdkPixbuf-2.0.typelib
exit 0
