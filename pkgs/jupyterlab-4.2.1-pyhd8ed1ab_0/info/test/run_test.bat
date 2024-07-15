



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter lab --version
IF %ERRORLEVEL% NEQ 0 exit /B 1
jlpm --version
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter labextension list
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter lab licenses
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter lab path
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter server extension list
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter server extension list 1>server_extensions 2>&1
IF %ERRORLEVEL% NEQ 0 exit /B 1
grep -iE "jupyterlab.*4\.2\.1.*OK" server_extensions
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter lab build --dev-build=False --minimize=False
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter lab clean
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
