

set -ex



pip check
coverage run --source jsonschema_specifications --branch -m pytest -vv --color=yes --tb=long --pyargs jsonschema_specifications
coverage report --show-missing --skip-covered --fail-under=96
exit 0
