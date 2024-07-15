



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
pybtex --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
pybtex-convert --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
pybtex-format --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
