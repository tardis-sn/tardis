:: Trailing semicolon in this variable as set by current (2017/01)
:: conda-build breaks us. Manual fix:
set "MSYS2_ARG_CONV_EXCL=/AI;/AL;/OUT;/out"
:: Delegate to the Unixy script. We need to translate the key path variables
:: to be Unix-y rather than Windows-y, though.
copy "%RECIPE_DIR%\build.sh" .
FOR /F "delims=" %%i IN ('cygpath.exe -u -p "%PATH%"') DO set "PATH_OVERRIDE=%%i"
FOR /F "delims=" %%i IN ('cygpath.exe -u "%LIBRARY_PREFIX%"') DO set "PREFIX=%%i"
FOR /F "delims=" %%i in ('cygpath.exe -u "%BUILD_PREFIX%"') DO set "BUILD_PREFIX=%%i"
set "SHLIB_EXT=.lib"
set "CMAKE_GENERATOR=MSYS Makefiles"
set MSYSTEM=MINGW%ARCH%
set MSYS2_PATH_TYPE=inherit
set CHERE_INVOKING=1
set "SHLIB_PREFIX="
set "fortran_compiler=m2w64-toolchain"
set "fortran_compiler_version=2"
bash -x "./build.sh"
IF %ERRORLEVEL% NEQ 0 exit 1
exit 0
