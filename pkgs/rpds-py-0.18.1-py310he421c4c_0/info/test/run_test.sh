

set -ex



pip check
coverage run --source=rpds -m pytest -vv --tb=long --color=yes
coverage report --show-missing --skip-covered --fail-under=100
exit 0
