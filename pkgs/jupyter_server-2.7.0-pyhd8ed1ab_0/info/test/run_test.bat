



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter server -h
IF %ERRORLEVEL% NEQ 0 exit /B 1
pytest -vv --cov=jupyter_server -k "not (delete or merge_config)" --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=70
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
