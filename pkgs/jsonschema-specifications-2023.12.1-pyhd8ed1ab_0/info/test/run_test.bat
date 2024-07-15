



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
coverage run --source jsonschema_specifications --branch -m pytest -vv --color=yes --tb=long --pyargs jsonschema_specifications
IF %ERRORLEVEL% NEQ 0 exit /B 1
coverage report --show-missing --skip-covered --fail-under=96
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
