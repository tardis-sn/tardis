echo on

mkdir build
cd build

REM remove GL flag for now
set "CXXFLAGS=-MD"
set "CC=cl.exe"
set "CXX=cl.exe"

cmake -G "Ninja" ^
    -DCMAKE_BUILD_TYPE="Release" ^
    -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
    -DCMAKE_INSTALL_PREFIX:PATH=%LIBRARY_PREFIX% ^
    -DLLVM_USE_INTEL_JITEVENTS=ON ^
    -DLLVM_ENABLE_LIBXML2=ON ^
    -DLLVM_ENABLE_RTTI=ON ^
    -DLLVM_ENABLE_ZLIB=ON ^
    -DLLVM_INCLUDE_BENCHMARKS=OFF ^
    -DLLVM_INCLUDE_DOCS=OFF ^
    -DLLVM_INCLUDE_EXAMPLES=OFF ^
    -DLLVM_INCLUDE_TESTS=OFF ^
    -DLLVM_INCLUDE_UTILS=ON ^
    -DLLVM_INSTALL_UTILS=ON ^
    -DLLVM_UTILS_INSTALL_DIR=libexec\llvm ^
    -DLLVM_BUILD_LLVM_C_DYLIB=no ^
    -DLLVM_EXPERIMENTAL_TARGETS_TO_BUILD=WebAssembly ^
    %SRC_DIR%/llvm
if %ERRORLEVEL% neq 0 exit 1

cmake --build .
if %ERRORLEVEL% neq 0 exit 1

:: SVML regression test (cf. #52 and numba/numba#3016)
bin\opt -S -vector-library=SVML -mcpu=haswell -O3 %RECIPE_DIR%\numba-3016.ll | bin\FileCheck %RECIPE_DIR%\numba-3016.ll
if %ERRORLEVEL% neq 0 exit 1

cd ..\llvm\test
..\..\build\bin\llvm-lit.py -vv Transforms ExecutionEngine Analysis CodeGen/X86
