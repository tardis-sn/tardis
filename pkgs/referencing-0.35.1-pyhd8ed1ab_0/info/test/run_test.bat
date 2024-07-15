



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
pytest -vv --pyargs referencing --cov=referencing --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=99
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
