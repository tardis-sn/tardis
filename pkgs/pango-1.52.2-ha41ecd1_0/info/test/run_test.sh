

set -ex



pango-view --help
test -f $PREFIX/lib/libpango-1.0${SHLIB_EXT}
test -f `pkg-config --variable=libdir --dont-define-prefix pango`/libpango-1.0${SHLIB_EXT}
test -f $PREFIX/lib/girepository-1.0/Pango-1.0.typelib
test -f $PREFIX/lib/libpangocairo-1.0${SHLIB_EXT}
test -f `pkg-config --variable=libdir --dont-define-prefix pangocairo`/libpangocairo-1.0${SHLIB_EXT}
test -f $PREFIX/lib/girepository-1.0/PangoCairo-1.0.typelib
test -f $PREFIX/lib/libpangoft2-1.0${SHLIB_EXT}
test -f `pkg-config --variable=libdir --dont-define-prefix pangoft2`/libpangoft2-1.0${SHLIB_EXT}
test -f $PREFIX/lib/girepository-1.0/PangoFT2-1.0.typelib
exit 0
