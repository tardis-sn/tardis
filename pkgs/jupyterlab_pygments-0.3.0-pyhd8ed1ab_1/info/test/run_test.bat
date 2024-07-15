



test -f ${PREFIX}/share/jupyter/labextensions/jupyterlab_pygments/package.json
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
