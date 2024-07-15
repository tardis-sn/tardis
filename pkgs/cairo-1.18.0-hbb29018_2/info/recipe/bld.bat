@ECHO ON

set "PKG_CONFIG_PATH=%LIBRARY_LIB%\pkgconfig;%LIBRARY_PREFIX%\share\pkgconfig;%BUILD_PREFIX%\Library\lib\pkgconfig"

:: get the prefix in "mixed" form
set "LIBRARY_PREFIX_M=%LIBRARY_PREFIX:\=/%"

%BUILD_PREFIX%\Scripts\meson setup builddir ^
  --buildtype=release ^
  --default-library=both ^
  --prefix=%LIBRARY_PREFIX_M% ^
  --wrap-mode=nofallback ^
  --backend=ninja ^
  -Dfontconfig=enabled ^
  -Dfreetype=enabled ^
  -Dglib=enabled
if errorlevel 1 exit 1

ninja -v -C builddir -j %CPU_COUNT%
if errorlevel 1 exit 1

ninja -C builddir install -j %CPU_COUNT%
if errorlevel 1 exit 1

move %LIBRARY_LIB%\libcairo.a %LIBRARY_LIB%\cairo-static.lib
move %LIBRARY_LIB%\libcairo-gobject.a %LIBRARY_LIB%\cairo-gobject-static.lib
move %LIBRARY_LIB%\libcairo-script-interpreter.a %LIBRARY_LIB%\cairo-script-interpreter-static.lib
