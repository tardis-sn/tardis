

set -ex



pip check
cd src/tests && pytest --cov=fqdn --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=95
exit 0
