

set -ex



test -f "${PREFIX}/include/yaml.h"
test -f "${PREFIX}/lib/libyaml${SHLIB_EXT}"
test -f "${PREFIX}/lib/pkgconfig/yaml-0.1.pc"
exit 0
