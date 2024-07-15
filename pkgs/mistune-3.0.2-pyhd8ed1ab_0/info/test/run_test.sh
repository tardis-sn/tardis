

set -ex



pip check
pytest -vv --color=yes --tb=long --cov=mistune --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=95
exit 0
