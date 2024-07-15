setlocal EnableDelayedExpansion
@echo on

set "PKG_CONFIG_PATH=%LIBRARY_LIB%\pkgconfig;%LIBRARY_PREFIX%\share\pkgconfig;%BUILD_PREFIX%\Library\lib\pkgconfig"

set ^"MESON_OPTIONS=^
  --prefix="%LIBRARY_PREFIX%" ^
  --default-library=shared ^
  --wrap-mode=nofallback ^
  --buildtype=release ^
  --backend=ninja ^
 ^"

meson setup builddir !MESON_OPTIONS!
if errorlevel 1 exit 1

meson configure builddir
if errorlevel 1 exit 1

ninja -v -C builddir -j %CPU_COUNT%
if errorlevel 1 exit 1

ninja -v -C builddir test -j %CPU_COUNT%
if errorlevel 1 exit 1

ninja -C builddir install -j %CPU_COUNT%
if errorlevel 1 exit 1
