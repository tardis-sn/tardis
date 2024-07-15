



test -f ${PREFIX}/fonts/SourceCodePro-Black.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-BlackIt.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-Bold.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-BoldIt.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-ExtraLight.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-ExtraLightIt.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-It.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-Light.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-LightIt.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-Medium.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-MediumIt.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-Regular.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-Semibold.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
test -f ${PREFIX}/fonts/SourceCodePro-SemiboldIt.ttf
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
