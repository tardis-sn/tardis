pushd "%SRC_DIR%"\build_cmake
ninja install
if errorlevel 1 exit 1

mkdir "%SRC_DIR%"\static_libs_for_cf
move /Y "%LIBRARY_LIB%"\pcre2-*-static.lib "%SRC_DIR%"\static_libs_for_cf\
if errorlevel 1 exit 1
