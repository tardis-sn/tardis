:: Need to set manifest for VS2008 otherwise exes crash (not sure why).
set MANIFEST=no
if "%VS_MAJOR%" == "9" (
    set MANIFEST=yes
)

:: Get libxml2 include paths via pkg-config
:: N.B.: We may at some point want to try another build system that supports
::       pkg-config out of the box.
set "PKG_CONFIG_PATH=%LIBRARY_LIB%\pkgconfig;%LIBRARY_PREFIX%\share\pkgconfig"
for /F "usebackq delims=" %%f in (`pkg-config --cflags-only-I libxml-2.0`) do set "libxml2_include=%%f"
set "libxml2_include=%libxml2_include: -I=;%"
set "libxml2_include=%libxml2_include:-I=;%"

cd win32
cscript configure.js prefix=%LIBRARY_PREFIX% include=%LIBRARY_INC%%libxml2_include% ^
        lib=%LIBRARY_LIB% sodir=%LIBRARY_BIN% iconv=yes zlib=yes vcmanifest=%MANIFEST%
if errorlevel 1 exit 1

nmake /f Makefile.msvc
if errorlevel 1 exit 1

nmake /f Makefile.msvc install
if errorlevel 1 exit 1
