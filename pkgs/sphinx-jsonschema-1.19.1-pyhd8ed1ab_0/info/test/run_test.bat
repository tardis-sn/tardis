



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
python -c "import importlib; importlib.import_module('sphinx-jsonschema')"
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
