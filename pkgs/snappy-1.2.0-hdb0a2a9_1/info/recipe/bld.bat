setlocal EnableDelayedExpansion

:: Remove -GL from CXXFLAGS as this causes a fatal error
set "CFLAGS= -MD"
set "CXXFLAGS= -MD"

mkdir build-dynamic
cd build-dynamic

cmake -G "NMake Makefiles" ^
    -DCMAKE_INSTALL_PREFIX:PATH="%LIBRARY_PREFIX%" ^
    -DCMAKE_PREFIX_PATH:PATH="%LIBRARY_PREFIX%" ^
    -DCMAKE_BUILD_TYPE:STRING=Release ^
    -DSNAPPY_BUILD_TESTS=OFF ^
    -DSNAPPY_BUILD_BENCHMARKS=OFF ^
    -DBUILD_SHARED_LIBS=ON ^
    ..
if errorlevel 1 exit 1

nmake
if errorlevel 1 exit 1

nmake install
if errorlevel 1 exit 1

cd ..
