setlocal EnableDelayedExpansion
@echo on

cd "win32"

:: set pkg-config path so that host deps can be found
set "PKG_CONFIG_PATH=%LIBRARY_LIB%\pkgconfig;%LIBRARY_PREFIX%\share\pkgconfig"

:: set XDG_DATA_DIRS to find gir files
set "XDG_DATA_DIRS=%XDG_DATA_DIRS%;%LIBRARY_PREFIX%\share"

:: add include dirs to search path
set "INCLUDE=%INCLUDE%;%LIBRARY_INC%\cairo;%LIBRARY_INC%\gdk-pixbuf-2.0"

:: build options
:: (override rustup command so that the conda-forge rust installation is used)
:: (add libiconv for linking against because glib needs its symbols)
:: (abuse LIBINTL_LIB to add libs that are needed for linking RSVG tools)
:: (override BINDIR to ensure the gobject-introspection tools are found)
set ^"LIBRSVG_OPTIONS=^
  CFG=release ^
  PREFIX="%LIBRARY_PREFIX%" ^
  BINDIR="%BUILD_PREFIX%\Library\bin" ^
  INTROSPECTION=1 ^
  RUSTUP=echo ^
  PYTHON="%BUILD_PREFIX%\python.exe" ^
  TOOLCHAIN_TYPE=stable ^
  LIBINTL_LIB="intl.lib iconv.lib advapi32.lib bcrypt.lib ws2_32.lib userenv.lib ntdll.lib" ^
  CARGO_CMD="cargo --locked build --release $(MANIFEST_PATH_FLAG) $(CARGO_TARGET_DIR_FLAG)" ^
 ^"

:: configure files
:: (use cmake just because it's convenient for replacing @VAR@ in files
cmake -DPACKAGE_VERSION=%PKG_VERSION% -P "%RECIPE_DIR%\win_configure_files.cmake"
if errorlevel 1 exit 1

nmake /F Makefile.vc !LIBRSVG_OPTIONS!
if errorlevel 1 exit 1

nmake /F Makefile.vc install !LIBRSVG_OPTIONS!
if errorlevel 1 exit 1

:: don't include debug symbols
del %LIBRARY_BIN%\rsvg-*.pdb
if errorlevel 1 exit 1
del %LIBRARY_LIB%\gdk-pixbuf-2.0\2.10.0\loaders\libpixbufloader-svg.pdb
if errorlevel 1 exit 1
