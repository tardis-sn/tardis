setlocal EnableDelayedExpansion
echo on

cmake -LAH -S .                               ^
   -B build_conda -G "Ninja"                  ^
   -DCMAKE_BUILD_TYPE="Release"               ^
   -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%"
if errorlevel 1 exit 1

cmake --build build_conda -- install
if errorlevel 1 exit 1
