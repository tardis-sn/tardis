@echo on
if "%ARCH%" == "64" (
   set ARCH=x64
) else (
   set ARCH=Win32
)

:: The vc9 build uses a patch that added CMake support from a fork.
:: We do not want to use it for vc14 since upstream provide a sln.

:: Let us not bother with the bundled VS solution files since they
:: do not build the CLI programs.
:: if "%vc%" == "9" goto vc9_build
:: goto vc14_build

cmake -GNinja ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
      -DCMAKE_C_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=FALSE ^
      -DCMAKE_BUILD_TYPE=Release ^
      --debug-trycompile ^
      .
ninja
if errorlevel 1 exit /b 1
