:: Make sure that gdk-pixbuf's loaders.cache is fully up-to-date.

:: Packages that install gdk-pixbuf loaders (such as librsvg) should have
:: post-link and post-unlink scripts that just execute this one, which will be
:: available as `%PREFIX%\Scripts\.gdk-pixbuf-post-link.bat`.

"%PREFIX%\Library\bin\gdk-pixbuf-query-loaders.exe" --update-cache >> "%PREFIX%/.messages.txt" 2>&1
if errorlevel 1 exit 1
