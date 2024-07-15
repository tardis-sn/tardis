



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter nbconvert --version
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter dejavu --version
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter nbconvert --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter dejavu --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
