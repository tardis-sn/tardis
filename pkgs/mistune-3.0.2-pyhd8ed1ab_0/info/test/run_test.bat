



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
pytest -vv --color=yes --tb=long --cov=mistune --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=95
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
