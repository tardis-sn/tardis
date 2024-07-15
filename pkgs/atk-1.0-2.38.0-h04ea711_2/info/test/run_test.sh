

set -ex



test -f $PREFIX/lib/libatk-1.0${SHLIB_EXT}
test -f `pkg-config --variable=libdir --dont-define-prefix atk`/libatk-1.0${SHLIB_EXT}
exit 0
