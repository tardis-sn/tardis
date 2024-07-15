setlocal EnableDelayedExpansion
rmdir /s /q internal-complibs

mkdir build
if errorlevel 1 exit 1
cd build
if errorlevel 1 exit 1

cmake -G "NMake Makefiles" ^
      %CMAKE_ARGS% ^
      -DCMAKE_BUILD_TYPE:STRING="Release" ^
      -DCMAKE_INSTALL_PREFIX:PATH="%LIBRARY_PREFIX%" ^
      -DCMAKE_PREFIX_PATH:PATH="%LIBRARY_PREFIX%" ^
      -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON ^
      -DBUILD_STATIC:BOOL=OFF ^
      -DBUILD_SHARED:BOOL=ON ^
      -DBUILD_TESTS:BOOL=ON ^
      -DBUILD_EXAMPLES:BOOL=OFF ^
      -DBUILD_BENCHMARKS:BOOL=OFF ^
      -DBLOSC_IS_SUBPROJECT=ON ^
      -DBLOSC_INSTALL=ON ^
      -DPREFER_EXTERNAL_LZ4:BOOL=ON ^
      -DPREFER_EXTERNAL_ZSTD:BOOL=ON ^
      -DPREFER_EXTERNAL_ZLIB:BOOL=ON ^
      "%SRC_DIR%"
if errorlevel 1 exit 1

cmake --build . --config Release
if errorlevel 1 exit 1

ctest -C release --output-on-failure --timeout 120
if errorlevel 1 exit 1

cmake --build . --target install --config Release
if errorlevel 1 exit 1
