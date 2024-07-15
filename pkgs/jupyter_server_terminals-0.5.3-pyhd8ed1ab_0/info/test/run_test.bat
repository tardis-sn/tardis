



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
coverage run --source=jupyter_server_terminals --branch -m pytest -vv --color=yes --tb=long
IF %ERRORLEVEL% NEQ 0 exit /B 1
coverage report --show-missing --skip-covered --fail-under=80
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
