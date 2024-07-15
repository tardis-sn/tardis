
:: Hack for winres.h being called winresrc.h on VS2008
if %VS_MAJOR% LEQ 9 copy %RECIPE_DIR%\winres.h .
if errorlevel 1 exit 1

nmake /f Makefile.vc
if errorlevel 1 exit 1

:: Build!
nmake /f Makefile.vc CFG=release-dynamic RTLIBCFG=dynamic OBJDIR=output UNICODE=0 all
if errorlevel 1 exit 1

:: gif2webp is build separately and not included in the "all" target
:: nmake /f Makefile.vc CFG=release-dynamic RTLIBCFG=dynamic OBJDIR=output UNICODE=0 gif2webp
:: if errorlevel 1 exit 1

copy output\release-dynamic\%ARCH%\bin\cwebp.exe %LIBRARY_PREFIX%\bin\cwebp.exe
if errorlevel 1 exit 1
copy output\release-dynamic\%ARCH%\bin\dwebp.exe %LIBRARY_PREFIX%\bin\dwebp.exe
if errorlevel 1 exit 1
:: copy output\release-dynamic\%ARCH%\bin\gif2webp.exe %LIBRARY_PREFIX%\bin\gif2webp.exe
:: if errorlevel 1 exit 1
copy output\release-dynamic\%ARCH%\bin\img2webp.exe %LIBRARY_PREFIX%\bin\img2webp.exe
if errorlevel 1 exit 1
copy output\release-dynamic\%ARCH%\bin\webpinfo.exe %LIBRARY_PREFIX%\bin\webpinfo.exe
if errorlevel 1 exit 1
copy output\release-dynamic\%ARCH%\bin\webpmux.exe %LIBRARY_PREFIX%\bin\webpmux.exe
if errorlevel 1 exit 1
