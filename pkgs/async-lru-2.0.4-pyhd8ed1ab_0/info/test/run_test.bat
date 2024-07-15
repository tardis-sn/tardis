



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
mypy -m async_lru
IF %ERRORLEVEL% NEQ 0 exit /B 1
pytest -vv --asyncio-mode=auto --cov=async_lru --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=89
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
