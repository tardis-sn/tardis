



pip check
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter kernelspec list
IF %ERRORLEVEL% NEQ 0 exit /B 1
jupyter run -h
IF %ERRORLEVEL% NEQ 0 exit /B 1
coverage run --source=jupyter_client --branch -m pytest -vv --tb=long --color=yes -k "not (signal_kernel_subprocesses or start_parallel_thread_kernels or start_parallel_process_kernels or open_tunnel or load_ips or input_request or tcp_cinfo)"
IF %ERRORLEVEL% NEQ 0 exit /B 1
coverage report --show-missing --skip-covered --fail-under=73
IF %ERRORLEVEL% NEQ 0 exit /B 1
exit /B 0
