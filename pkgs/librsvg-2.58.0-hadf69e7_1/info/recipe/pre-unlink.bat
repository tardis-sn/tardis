set GDK_PIXBUF_POST_LINK_SCRIPT=%PREFIX%\Scripts\.gdk-pixbuf-post-link.bat
if not exist "%GDK_PIXBUF_POST_LINK_SCRIPT%" exit 0

:: Since librsvg is being removed, we want it to be removed from loaders.cache.
:: But since this is a PRE-unlink script, the loader is still present.
:: Remove it now so we can update loaders.cache correctly.
cd "%PREFIX%\Library\lib\gdk-pixbuf-2.0" > nul 2>> "%PREFIX%/.messages.txt"
if errorlevel 1 exit 1

del /F /S /Q "libpixbufloader-svg.dll" > nul 2>> "%PREFIX%/.messages.txt"
if errorlevel 1 exit 1

:: The gdk-pixbuf post-link function updates the loaders
call "%GDK_PIXBUF_POST_LINK_SCRIPT%"
