



pytest -vv --cov=anyio --cov-branch --cov-report=term-missing:skip-covered -k="not (bind_link_local or block_device or connection_refused or happy_eyeballs or ipv6)" --cov-fail-under=86
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
