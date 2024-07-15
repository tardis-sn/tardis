:: Primary build
set _gz_builddir=%SRC_DIR%\build
cmake -S%SRC_DIR% ^
      -GNinja ^
      -B%_gz_builddir% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
      -DCMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP=ON ^
      -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX%;%LIBRARY_LIB% ^
      -DLTDL_INCLUDE_DIR=%SRC_DIR%\ltdl_compat ^
      -Dwith_gvedit=OFF ^
      -Duse_win_pre_inst_libs=OFF ^
      -Dinstall_win_dependency_dlls=OFF
if errorlevel 1 exit 1

cmake --build build -- install
if errorlevel 1 exit 1


:: Reinstall into a temporary directory to gather names of execs
set _gz_installdir=%SRC_DIR%\install
cmake -Wno-dev -DCMAKE_INSTALL_PREFIX=%_gz_installdir% %_gz_builddir%
if errorlevel 1 exit 1

cmake --build %_gz_builddir% -- install
if errorlevel 1 exit 1


:: Setup wrappers for backwards compatibility
if not exist "%PREFIX%\Scripts" mkdir %PREFIX%\Scripts
cd %PREFIX%\Scripts
for /r "%_gz_installdir%\bin" %%f in (*.exe) do (
    echo @echo off > %%~nf.bat
    echo "%%~dp0.\..\Library\bin\%%~nf.exe" %%* >> %%~nf.bat
    if errorlevel 1 exit 1
)
