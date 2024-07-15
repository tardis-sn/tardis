REM Build the lib...
nmake /F %RECIPE_DIR%\Makefile.vc build_libs
if errorlevel 1 exit 1

REM done compiling

REM Install step
copy src\*.h %LIBRARY_INC%\	
if errorlevel 1 exit 1
cd .\gdbuild
copy libgd.dll %LIBRARY_BIN%\
if errorlevel 1 exit 1
copy libgd.lib %LIBRARY_LIB%\
if errorlevel 1 exit 1
copy libgd.lib %LIBRARY_LIB%\gd.lib
if errorlevel 1 exit 1
copy libgd_a.lib %LIBRARY_LIB%\
if errorlevel 1 exit 1