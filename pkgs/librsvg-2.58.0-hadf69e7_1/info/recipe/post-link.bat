set GDK_PIXBUF_POST_LINK_SCRIPT=%PREFIX%\Scripts\.gdk-pixbuf-post-link.bat

:: The gdk-pixbuf post-link function updates the loaders
if exist "%GDK_PIXBUF_POST_LINK_SCRIPT%" (
    call "%GDK_PIXBUF_POST_LINK_SCRIPT%"
) else >> "%PREFIX%\.messages.txt" (
    echo.librsvg: The post-link script did not complete.
    echo.To take advantage of gdk-pixbuf's support for librsvg, please run:
    echo.    %GDK_PIXBUF_POST_LINK_SCRIPT%
)
