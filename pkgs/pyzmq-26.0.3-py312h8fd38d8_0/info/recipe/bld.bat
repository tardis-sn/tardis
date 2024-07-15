set PYZMQ_NO_BUNDLE=1

"%PYTHON%" -m pip install .
if errorlevel 1 exit 1
