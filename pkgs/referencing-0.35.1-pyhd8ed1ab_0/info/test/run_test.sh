

set -ex



pip check
pytest -vv --pyargs referencing --cov=referencing --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=99
exit 0
