



test -f ${PREFIX}/fonts/Ubuntu-B.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-BI.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-C.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-L.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-LI.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-M.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-MI.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-R.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-RI.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/Ubuntu-Th.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/UbuntuMono-B.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/UbuntuMono-BI.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/UbuntuMono-R.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/UbuntuMono-RI.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
