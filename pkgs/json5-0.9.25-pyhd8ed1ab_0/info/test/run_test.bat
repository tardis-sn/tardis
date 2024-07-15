



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
pyjson5 --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
