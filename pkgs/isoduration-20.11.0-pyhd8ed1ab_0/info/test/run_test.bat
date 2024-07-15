



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
cd src/tests && pytest --cov=isoduration --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=98
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
