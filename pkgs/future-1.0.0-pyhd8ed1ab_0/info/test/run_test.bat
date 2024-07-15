



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
futurize --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
pasteurize --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
