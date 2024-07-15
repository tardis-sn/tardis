

set -ex



pip check
jupyter kernelspec list
jupyter run -h
coverage run --source=jupyter_client --branch -m pytest -vv --tb=long --color=yes -k "not (signal_kernel_subprocesses or start_parallel_thread_kernels or start_parallel_process_kernels or open_tunnel or load_ips or input_request or tcp_cinfo)"
coverage report --show-missing --skip-covered --fail-under=73
exit 0
