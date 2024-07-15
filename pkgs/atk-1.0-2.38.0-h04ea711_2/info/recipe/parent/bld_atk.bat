setlocal EnableDelayedExpansion
@echo on

:: set pkg-config path so that host deps can be found
:: (set as env var so it's used by both meson and during build with g-ir-scanner)
set "PKG_CONFIG_PATH=%LIBRARY_LIB%\pkgconfig;%LIBRARY_PREFIX%\share\pkgconfig;%BUILD_PREFIX%\Library\lib\pkgconfig"

:: get mixed path (forward slash) form of prefix so host prefix replacement works
set "LIBRARY_PREFIX_M=%LIBRARY_PREFIX:\=/%"

:: meson options
:: (set pkg_config_path so deps in host env can be found)
:: (link to intl.lib since glib doesn't anymore and we need it)
set ^"MESON_OPTIONS=^
  --prefix="%LIBRARY_PREFIX_M%" ^
  --wrap-mode=nofallback ^
  --buildtype=release ^
  --backend=ninja ^
  -D docs=false ^
  -D c_link_args=intl.dll.lib ^
 ^"

:: configure build using meson
meson setup builddir !MESON_OPTIONS!
if errorlevel 1 exit 1

:: print results of build configuration
meson configure builddir
if errorlevel 1 exit 1

ninja -v -C builddir -j %CPU_COUNT%
if errorlevel 1 exit 1

ninja -C builddir install -j %CPU_COUNT%
if errorlevel 1 exit 1
