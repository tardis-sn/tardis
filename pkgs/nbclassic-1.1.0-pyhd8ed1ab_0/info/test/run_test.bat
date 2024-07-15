



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter nbclassic --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter nbclassic-extension --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter nbclassic-serverextension --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter nbclassic-bundlerextension --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter server extension list
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter server extension list 1>server_extensions 2>&1
IF %ERRORLEVEL% NEQ 0 exit /B 1
cat server_extensions | grep -iE "nbclassic.*enabled"
IF %ERRORLEVEL% NEQ 0 exit /B 1
cat server_extensions | grep -iE "nbclassic.*OK"
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
