



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
pytest --cov=rfc3339_validator --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=100
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
