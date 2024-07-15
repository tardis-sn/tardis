# Taken from https://sourceforge.net/p/giflib/feature-requests/6/
# Original Author: Sandro Mani

cmake_minimum_required(VERSION 2.6.0)

project(giflib C)

SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
SET(BUILD_STATIC_LIBS OFF CACHE BOOL "Whether to also build static libs")

execute_process(COMMAND ./getversion
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    message(FATAL_ERROR "In tree builds are not allowed. Use a build directory.")
endif()

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/Makefile" MAKEFILE_CONTENTS)
string(REGEX MATCH "LIBMAJOR=([0-9]*)" _ ${MAKEFILE_CONTENTS})
set(LIBMAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "LIBMINOR=([0-9]*)" _ ${MAKEFILE_CONTENTS})
set(LIBMINOR ${CMAKE_MATCH_1})
string(REGEX MATCH "LIBPOINT=([0-9]*)" _ ${MAKEFILE_CONTENTS})
set(LIBPOINT ${CMAKE_MATCH_1})

set(LIBVER "${LIBMAJOR}.${LIBMINOR}.${LIBPOINT}")

set(giflib_SRC
    dgif_lib.c
    egif_lib.c
    getarg.c
    gifalloc.c
    gif_err.c
    gif_font.c
    gif_hash.c
    openbsd-reallocarray.c
    qprintf.c
    quantize.c
)

# Some utilities are installed
set(giflib_INSTALLABLE
	gif2rgb
	gifbuild
	giffix
	giftext
	giftool
	gifclrmp
)

# Some utilities are only used internally for testing.
# There is a parallel list in doc/Makefile.
# These are all candidates for removal in future releases.
set(giflib_UTILS
    ${giflib_INSTALLABLE}
	gifbg
	gifcolor
	gifecho
	giffilter
	gifhisto
	gifinto
	gifsponge
	gifwedge
)

file(GLOB giflib_MAN doc/*.1)

### Build library / tools

add_library(gif SHARED ${giflib_SRC})

if (NOT WIN32)
    target_link_libraries(gif m)
endif()

set_target_properties(gif PROPERTIES VERSION ${LIBVER} SOVERSION ${LIBMAJOR})
if(WIN32)
    set_target_properties(gif
        PROPERTIES
        RUNTIME_OUTPUT_NAME "gif-${LIBMAJOR}"
        COMPILE_DEFINITIONS "GIFLIB_EXPORTS"
    )
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

if(${BUILD_STATIC_LIBS})
    add_library(gif_static STATIC ${giflib_SRC})
    if(WIN32)
        set_target_properties(gif_static PROPERTIES OUTPUT_NAME libgif)
    else()
        set_target_properties(gif_static PROPERTIES OUTPUT_NAME gif)
    endif()
endif(${BUILD_STATIC_LIBS})

### Installation

install(TARGETS gif
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib${LIB_SUFFIX}
    LIBRARY DESTINATION lib${LIB_SUFFIX}
)

if(${BUILD_STATIC_LIBS})
    install(TARGETS gif_static ARCHIVE DESTINATION lib${LIB_SUFFIX})
endif(${BUILD_STATIC_LIBS})

foreach(UTILITY ${giflib_UTILS})
    add_executable(${UTILITY} ${UTILITY}.c)
    target_link_libraries(${UTILITY} gif)
endforeach()

foreach(UTILITY ${giflib_INSTALLABLE})
    install(TARGETS ${UTILITY} DESTINATION bin)
endforeach()

install(FILES gif_lib.h DESTINATION include)
install(FILES ${giflib_MAN} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/man/man1)

### Distribution tarball
set(giflib_DIST
    *.c
    *.h
    README
    NEWS
    TODO
    COPYING
    getversion
    ChangeLog
    CMakeLists.txt
    build.adoc
    history.adoc
    control
    doc/whatsinagif
    doc/*.1
    doc/*.xml
    doc/*.txt
    doc/index.html.in
    doc/00README
    doc/Makefile
    tests
    pic
)
# We include all of the XML, and also generated manual pages
# so people working from the distribution tarball won't need xmlto.
add_custom_target(dist-gz
    COMMAND tar --transform='s:^:giflib-${VERSION}/:' -czf giflib-${VERSION}.tar.gz ${giflib_DIST}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_custom_target(dist-bz2
    COMMAND tar --transform='s:^:giflib-${VERSION}/:' -cjf giflib-${VERSION}.tar.bz2 ${giflib_DIST}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_custom_target(dist DEPENDS dist-gz DEPENDS dist-bz2)


### Auditing tools

# cppcheck should run clean
add_custom_target(cppcheck
    COMMAND cppcheck --inline-suppr --template gcc --enable=all --suppress=unusedFunction --force *.[ch]
)
