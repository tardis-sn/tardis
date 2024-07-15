@ECHO ON

:: set pkg-config path so that host deps can be found
:: (set as env var so it's used by both meson and during build with g-ir-scanner)
set "PKG_CONFIG_PATH=%LIBRARY_LIB%\pkgconfig;%LIBRARY_PREFIX%\share\pkgconfig;%BUILD_PREFIX%\Library\lib\pkgconfig"

:: get mixed path (forward slash) form of prefix so host prefix replacement works
set "LIBRARY_PREFIX_M=%LIBRARY_PREFIX:\=/%"

%BUILD_PREFIX%\Scripts\meson setup builddir ^
	--wrap-mode=nofallback ^
	--buildtype=release ^
	--prefix=%LIBRARY_PREFIX_M% ^
	--backend=ninja ^
	-Dglib=enabled ^
	-Dgobject=enabled ^
	-Dcairo=enabled ^
	-Dicu=enabled ^
	-Dgraphite2=enabled ^
	-Dfreetype=enabled ^
	-Dgdi=auto ^
	-Dcoretext=auto ^
	-Dintrospection=enabled ^
	-Dtests=disabled ^
	-Dbenchmark=disabled ^
	-Ddocs=disabled
if errorlevel 1 exit 1

ninja -v -C builddir -j %CPU_COUNT%
if errorlevel 1 exit 1

ninja -C builddir install -j %CPU_COUNT%
if errorlevel 1 exit 1

del %LIBRARY_PREFIX%\bin\*.pdb
