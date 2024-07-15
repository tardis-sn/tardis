

set -ex



pip check
cd src/tests && pytest --cov=isoduration --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=98
exit 0
