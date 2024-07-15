

set -ex



pip check
jupyter server -h
pytest -vv --cov=jupyter_server -k "not (delete or merge_config)" --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=70
exit 0
