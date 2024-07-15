@echo on

cd build
cmake --build . --target install
if %ERRORLEVEL% neq 0 exit 1

if "%PKG_NAME%" == "llvm-tools" (
    del /f %LIBRARY_BIN%\*.dll
    rmdir /S /Q %LIBRARY_LIB%
    rmdir /S /Q %LIBRARY_INC%
    rmdir /S /Q %LIBRARY_PREFIX%\libexec
) else (
    REM upstream picks up diaguids.lib from the windows image, see
    REM https://github.com/llvm/llvm-project/blob/llvmorg-14.0.6/llvm/lib/DebugInfo/PDB/CMakeLists.txt#L17
    REM which ultimately derives from VSINSTALLDIR, see
    REM https://github.com/llvm/llvm-project/blob/llvmorg-14.0.6/llvm/cmake/config-ix.cmake#L516
    REM and gets hardcoded by CMake to point to the path in our windows image.
    REM This makes it non-portable between image versions (e.g. 2019 vs 2022), so replace
    REM the hardcoded path with a variable again
    sed -i "s,C:/Program Files/Microsoft Visual Studio/2022/Enterprise,$ENV{VSINSTALLDIR},g" %LIBRARY_LIB%\cmake\llvm\LLVMExports.cmake
)
