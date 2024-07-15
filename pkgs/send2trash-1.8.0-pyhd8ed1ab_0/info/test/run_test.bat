



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
send2trash -h
IF %ERRORLEVEL% NEQ 0 exit /B 1
echo test > test.txt
IF %ERRORLEVEL% NEQ 0 exit /B 1
python -c "from send2trash import *; send2trash('test.txt')"
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
