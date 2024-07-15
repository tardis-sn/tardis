



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
pip list | grep -iE nbformat
IF %ERRORLEVEL% NEQ 0 exit /B 1
pip list | grep -iE 'nbformat\s*5\.10\.4'
IF %ERRORLEVEL% NEQ 0 exit /B 1
python -c "v = __import__('nbformat').__version__; print(v); assert v == '5.10.4'"
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter trust --version
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter-trust --help
IF %ERRORLEVEL% NEQ 0 exit /B 1
pytest -vv --cov=nbformat --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=76
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
