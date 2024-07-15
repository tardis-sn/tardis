



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
conda install -y numpy
IF %ERRORLEVEL% NEQ 0 exit /B 1
python test-imports-requiring-numpy.py
IF %ERRORLEVEL% NEQ 0 exit /B 1
conda install -y matplotlib-base
IF %ERRORLEVEL% NEQ 0 exit /B 1
python test-imports-requiring-matplotlib.py
IF %ERRORLEVEL% NEQ 0 exit /B 1
conda install -y pandas
IF %ERRORLEVEL% NEQ 0 exit /B 1
python test-imports-requiring-pandas.py
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
