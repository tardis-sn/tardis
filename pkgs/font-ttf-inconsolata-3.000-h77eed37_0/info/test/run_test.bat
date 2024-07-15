



test -f ${PREFIX}/fonts/Inconsolata-Regular.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Inconsolata-Bold.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
