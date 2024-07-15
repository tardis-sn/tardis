



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
cm2html --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
cm2latex --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
cm2man --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
cm2pseudoxml --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
cm2xetex --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
cm2xml --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
