# cmake script to do variable substitution for Windows build files
# PACKAGE_VERSION should be set when calling this
if(NOT DEFINED PACKAGE_VERSION)
    message(FATAL_ERROR "Must define PACKAGE_VERSION variable")
endif()

set(PACKAGE "RSVG")
set(PACKAGE_BUGREPORT "https://gitlab.gnome.org/GNOME/librsvg/issues")
set(PACKAGE_NAME ${PACKAGE})
set(PACKAGE_TARNAME "librsvg")
set(GETTEXT_PACKAGE "librsvg")

# extract major, minor, and micro from version string
string(REPLACE "." ";" VERSION_LIST ${PACKAGE_VERSION})
list(GET VERSION_LIST 0 LIBRSVG_MAJOR_VERSION)
list(GET VERSION_LIST 1 LIBRSVG_MINOR_VERSION)
list(GET VERSION_LIST 2 LIBRSVG_MICRO_VERSION)

configure_file(config-msvc.mak.in config-msvc.mak @ONLY)
configure_file(config.h.win32.in config.h.win32 @ONLY)

