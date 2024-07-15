

set -ex



pip check
coverage run --source=jupyter_server_terminals --branch -m pytest -vv --color=yes --tb=long
coverage report --show-missing --skip-covered --fail-under=80
exit 0
