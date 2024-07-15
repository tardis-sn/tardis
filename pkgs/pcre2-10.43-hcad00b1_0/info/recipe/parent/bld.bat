CFLAGS="${CFLAGS} -O3"
CXXFLAGS="${CXXFLAGS} -O3"

mkdir "%SRC_DIR%"\build_cmake
pushd "%SRC_DIR%"\build_cmake
cmake ^
    -DBUILD_SHARED_LIBS=ON ^
    -DBUILD_STATIC_LIBS=ON ^
    -DCMAKE_BUILD_TYPE=release ^
    -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
    -DPCRE2_SUPPORT_JIT=ON ^
    -DPCRE2_BUILD_PCRE2_16=ON ^
    -DPCRE2_BUILD_PCRE2_32=ON ^
    -GNinja ^
    ..
if errorlevel 1 exit 1

ninja
if errorlevel 1 exit 1
sed -ie "s/\$<TARGET_FILE:pcre2test>/pcre2test.exe/" pcre2_test.bat
cp -r ..\testdata .
if errorlevel 1 exit 1
ninja test
if errorlevel 1 exit 1
