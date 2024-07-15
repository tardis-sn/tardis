

set -ex



pcre2test --version
pcre2grep --version
pcre2-config --version
test -f ${PREFIX}/include/pcre2.h
test -f ${PREFIX}/include/pcre2posix.h
test -f ${PREFIX}/lib/libpcre2-posix${SHLIB_EXT}
test -f ${PREFIX}/lib/libpcre2-posix${SHLIB_EXT}.3
test -f ${PREFIX}/lib/libpcre2-posix${SHLIB_EXT}.3.0.5
test -f ${PREFIX}/lib/libpcre2-8${SHLIB_EXT}
test -f ${PREFIX}/lib/libpcre2-8${SHLIB_EXT}.0
test -f ${PREFIX}/lib/libpcre2-8${SHLIB_EXT}.0.13.0
test -f ${PREFIX}/lib/libpcre2-16${SHLIB_EXT}
test -f ${PREFIX}/lib/libpcre2-16${SHLIB_EXT}.0
test -f ${PREFIX}/lib/libpcre2-16${SHLIB_EXT}.0.13.0
test -f ${PREFIX}/lib/libpcre2-32${SHLIB_EXT}
test -f ${PREFIX}/lib/libpcre2-32${SHLIB_EXT}.0
test -f ${PREFIX}/lib/libpcre2-32${SHLIB_EXT}.0.13.0
test -f ${PREFIX}/lib/pkgconfig/libpcre2-8.pc
test -f ${PREFIX}/lib/pkgconfig/libpcre2-16.pc
test -f ${PREFIX}/lib/pkgconfig/libpcre2-32.pc
test -f ${PREFIX}/lib/pkgconfig/libpcre2-posix.pc
test -f ${PREFIX}/cmake/pcre2-config.cmake
test -f ${PREFIX}/cmake/pcre2-config-version.cmake
exit 0
