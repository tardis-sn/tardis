

set -ex



pip check
jupyter -h
jupyter-migrate -h
jupyter-troubleshoot --help
pytest -vv --color=yes --cov=jupyter_core --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=64 -k "not (test_not_on_path or test_path_priority or test_argv0)"
exit 0
