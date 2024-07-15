



python -c "from dataclasses import dataclass"
IF %ERRORLEVEL% NEQ 0 exit /B 1
python -c "import pkg_resources as p; p.get_distribution('dataclasses')"
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
