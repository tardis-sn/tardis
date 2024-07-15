

set -ex



test -f ${PREFIX}/lib/libpng.a
test -f ${PREFIX}/lib/libpng${SHLIB_EXT}
libpng-config --version
exit 0
