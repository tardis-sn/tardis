



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
pytest --cov=arrow --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=99 -k "not parse_tz_name_zzz"
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
