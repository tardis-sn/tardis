# Stop on first error
set -euxo pipefail

# test dir populated from recipe/test
pushd test

# Test C compiler
if [[ "$mpi" != "nompi" ]]; then
    h5cc=h5pcc
    h5fc=h5pfc
else
    h5cc=h5cc
    h5fc=h5fc
fi
echo "Testing $h5cc"
$h5cc -showconfig
$h5cc ${CFLAGS} ${LDFLAGS} h5_cmprss.c -o h5_cmprss
./h5_cmprss

# Test C++ compiler
echo "Testing h5c++"
h5c++ -showconfig
h5c++ ${CXXFLAGS} ${LDFLAGS} h5tutr_cmprss.cpp -o h5tutr_cmprss
./h5tutr_cmprss

# Test Fortran 90 compiler
echo "Testing $h5fc"
$h5fc -showconfig
$h5fc ${FFLAGS} ${LDFLAGS} h5_cmprss.f90 -o h5_cmprss
./h5_cmprss

# Test Fortran 2003 compiler, note that the file has a 90 extension
echo "Testing h5fc for Fortran 2003"
$h5fc ${FFLAGS} ${LDFLAGS} compound_fortran2003.f90 -o compound_fortran2003
./compound_fortran2003

echo "Testing cmake"
cmake -B build .
cmake --build build
./build/h5_cmprss
