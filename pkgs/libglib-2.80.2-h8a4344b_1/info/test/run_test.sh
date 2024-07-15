

set -ex



test -f ${PREFIX}/lib/libglib-2.0.so.0
test ! -f ${PREFIX}/lib/libgobject-2.0.la
test ! -f ${PREFIX}/lib/libglib-2.0${SHLIB_EXT}
test -f ${PREFIX}/lib/pkgconfig/glib-2.0.pc
test -f ${PREFIX}/etc/conda/activate.d/libglib_activate.sh
test -f ${PREFIX}/etc/conda/deactivate.d/libglib_deactivate.sh
test -f ${PREFIX}/share/gir-1.0/GLib-2.0.gir
test -f ${PREFIX}/lib/girepository-1.0/GLib-2.0.typelib
exit 0
