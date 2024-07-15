



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
python -c "import alabaster; print(alabaster.get_path())"
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
