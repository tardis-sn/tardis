mkdir build_%c_compiler%
cd build_%c_compiler%

cmake -GNinja                                    ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%"  ^
      -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%"     ^
      -DCMAKE_BUILD_TYPE=Release                 ^
      -DCMAKE_C_FLAGS="%CFLAGS% -DWIN32"         ^
      -DCMAKE_CXX_FLAGS="%CXXFLAGS% -EHsc"       ^
      -DCMAKE_SHARED_LIBRARY_PREFIX=""           ^
      ..
if errorlevel 1 exit /b 1

cmake --build . --target install --config Release
if errorlevel 1 exit /b 1

copy "%LIBRARY_PREFIX%"\bin\tiff.dll "%LIBRARY_PREFIX%"\bin\libtiff.dll
if errorlevel 1 exit /b 1
copy "%LIBRARY_PREFIX%"\lib\tiff.lib "%LIBRARY_PREFIX%"\lib\libtiff.lib
if errorlevel 1 exit /b 1

:REM https://gitlab.com/libtiff/libtiff/-/merge_requests/338
:REM copy "%LIBRARY_PREFIX%"\bin\tiffxx.dll "%LIBRARY_PREFIX%"\bin\libtiffxx.dll
:REM if errorlevel 1 exit /b 1
