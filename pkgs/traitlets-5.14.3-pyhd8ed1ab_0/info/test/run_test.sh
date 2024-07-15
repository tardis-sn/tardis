

set -ex



pip check
coverage run --source=traitlets --parallel --branch -m pytest -vv
coverage run --source=traitlets --parallel --branch -m pytest -vv --pyargs traitlets
coverage combine
coverage report --show-missing --skip-covered --fail-under=80
mypy -p traitlets
exit 0
