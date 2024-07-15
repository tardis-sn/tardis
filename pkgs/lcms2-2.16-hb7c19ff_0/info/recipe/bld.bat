set UseEnv=true

msbuild ^
  /p:Platform=x64 ^
  /p:Configuration=Release ^
  /p:AdditionalIncludeDirectories=%LIBRARY_INC% ^
  /p:AdditionalDependencies=/LIBPATH:%LIBRARY_LIB% ^
  Projects\VC2019\lcms2.sln
if errorlevel 1 exit 1

REM For debugging Purposes, you may want to list the files in the 3 important directories
REM dir
REM dir bin
REM dir include
REM dir lib

REM Note we renamed jpegicc and tifficc to match their linux counterparts
REM That existed in conda-forge before the windows versions.

COPY bin\jpegicc.exe  %LIBRARY_BIN%\jpgicc.exe
COPY bin\tifficc.exe  %LIBRARY_BIN%\tificc.exe
COPY bin\linkicc.exe  %LIBRARY_BIN%\linkicc.exe
COPY bin\transicc.exe %LIBRARY_BIN%\transicc.exe
COPY bin\psicc.exe    %LIBRARY_BIN%\psicc.exe

COPY bin\lcms2.dll    %LIBRARY_BIN%\lcms2.dll

COPY bin\lcms2.lib    %LIBRARY_LIB%\lcms2.lib

COPY include\lcms2.h         %LIBRARY_INC%\lcms2.h
COPY include\lcms2_plugin.h  %LIBRARY_INC%\lcms2_plugin.h
