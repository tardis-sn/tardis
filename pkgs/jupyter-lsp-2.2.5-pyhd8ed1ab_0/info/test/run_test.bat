



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter server extension list
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter server extension list 1>server_extensions 2>&1
IF %ERRORLEVEL% NEQ 0 exit /B 1
grep -iE "jupyter_lsp.*OK" server_extensions
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
