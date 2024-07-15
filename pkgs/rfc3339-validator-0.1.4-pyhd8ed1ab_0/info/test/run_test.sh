

set -ex



pip check
pytest --cov=rfc3339_validator --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=100
exit 0
