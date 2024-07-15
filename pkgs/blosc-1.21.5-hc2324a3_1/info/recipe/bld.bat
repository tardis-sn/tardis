setlocal EnableDelayedExpansion
rmdir /s /q internal-complibs

mkdir build
if errorlevel 1 exit 1
cd build
if errorlevel 1 exit 1

cmake -LAH -G "NMake Makefiles" ^
      -DCMAKE_BUILD_TYPE:STRING="Release" ^
      -DCMAKE_PREFIX_PATH:PATH="%LIBRARY_PREFIX%" ^
      -DCMAKE_INSTALL_PREFIX:PATH="%LIBRARY_PREFIX%" ^
      -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON ^
      -DBUILD_STATIC:BOOL=OFF ^
      -DBUILD_SHARED:BOOL=ON ^
      -DBUILD_TESTS:BOOL=ON ^
      -DBUILD_BENCHMARKS:BOOL=OFF ^
      -DDEACTIVATE_SNAPPY:BOOL=OFF ^
      -DPREFER_EXTERNAL_LZ4:BOOL=ON ^
      -DPREFER_EXTERNAL_ZLIB:BOOL=ON ^
      -DPREFER_EXTERNAL_ZSTD:BOOL=ON ^
      "%SRC_DIR%"
if errorlevel 1 exit 1

cmake --build . --config Release
if errorlevel 1 exit 1

ctest -C release
if errorlevel 1 exit 1

cmake --build . --target install --config Release
if errorlevel 1 exit 1

del %LIBRARY_BIN%\msvc*.dll
