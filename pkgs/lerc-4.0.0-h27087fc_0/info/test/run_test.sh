

set -ex



test -f ${PREFIX}/lib/libLerc${SHLIB_EXT}
test -f ${PREFIX}/include/Lerc_types.h
exit 0
