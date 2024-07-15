



pytest -v tests/ -k "not test_nested_typevar_in_signature"
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
