@ECHO ON

mkdir build
cd build

cmake -G Ninja ^
      -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DBUILD_SHARED_LIBS=ON ^
      -DGRAPHITE2_COMPARE_RENDERER=ON ^
      ..
if errorlevel 1 exit 1

ninja -v
if errorlevel 1 exit 1

ninja install
if errorlevel 1 exit 1
