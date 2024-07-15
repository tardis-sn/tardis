



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
pip list
IF %ERRORLEVEL% NEQ 0 exit /B 1
pip list | grep -iE "uri-template\s+1\.3\.0"
IF %ERRORLEVEL% NEQ 0 exit /B 1
mypy -p uri_template
IF %ERRORLEVEL% NEQ 0 exit /B 1
cd src && coverage run --source=uri_template test.py && coverage report -m --fail-under=91
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
