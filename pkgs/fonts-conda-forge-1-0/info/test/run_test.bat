



ls -alh ${PREFIX}/fonts/
IF %ERRORLEVEL% NEQ 0 exit /B 1
du -sh ${PREFIX}/fonts/
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
