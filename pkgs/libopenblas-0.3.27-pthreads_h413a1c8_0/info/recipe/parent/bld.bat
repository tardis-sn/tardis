mkdir build
cd build

cmake -G "NMake Makefiles JOM"              ^
    -DCMAKE_C_COMPILER=clang-cl             ^
    -DCMAKE_Fortran_COMPILER=flang          ^
    -DCMAKE_BUILD_TYPE=Release              ^
    -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
    -DDYNAMIC_ARCH=ON                       ^
    -DBUILD_WITHOUT_LAPACK=no               ^
    -DNOFORTRAN=0                           ^
    -DNUM_THREADS=128                       ^
    -DBUILD_SHARED_LIBS=on                  ^
    %SRC_DIR%

jom install -j%CPU_COUNT%

utest\openblas_utest.exe
