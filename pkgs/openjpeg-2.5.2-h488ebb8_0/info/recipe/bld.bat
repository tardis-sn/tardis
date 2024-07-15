mkdir build
cd build

cmake -GNinja ^
      -D CMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
      -D TIFF_LIBRARY=%LIBRARY_LIB%\tiff.lib ^
      -D TIFF_INCLUDE_DIR=%LIBRARY_INC% ^
      -D PNG_LIBRARY_RELEASE=%LIBRARY_LIB%\libpng.lib ^
      -D PNG_PNG_INCLUDE_DIR=%LIBRARY_INC% ^
      -D ZLIB_LIBRARY=%LIBRARY_LIB%\zlib.lib ^
      -D ZLIB_INCLUDE_DIR=%LIBRARY_INC% ^
      -D CMAKE_BUILD_TYPE=Release ^
      %SRC_DIR%
if errorlevel 1 exit 1

:: Build.
cmake --build . --config Release
if errorlevel 1 exit 1

:: Install.
cmake --build . --config Release --target install
if errorlevel 1 exit 1

:: Test.
ctest -C Release
if errorlevel 1 exit 1
