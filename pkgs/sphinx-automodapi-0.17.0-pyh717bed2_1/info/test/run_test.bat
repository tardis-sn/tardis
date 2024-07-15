



python -m pip check sphinx-automodapi
IF %ERRORLEVEL% NEQ 0 exit /B 1
python -m pip show sphinx-automodapi
IF %ERRORLEVEL% NEQ 0 exit /B 1
python -m pytest -ra --pyargs sphinx_automodapi -k 'not cython'
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
