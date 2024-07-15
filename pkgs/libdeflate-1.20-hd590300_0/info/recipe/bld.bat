@echo on

mkdir build
if errorlevel 1 exit 1

cd build
if errorlevel 1 exit 1


cmake -LAH -G "Ninja"                     ^
  %CMAKE_ARGS%                            ^
  -DLIBDEFLATE_BUILD_STATIC_LIB=OFF       ^
  -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
  ..
if errorlevel 1 exit 1

cmake --build . --target install --config Release
if errorlevel 1 exit 1
