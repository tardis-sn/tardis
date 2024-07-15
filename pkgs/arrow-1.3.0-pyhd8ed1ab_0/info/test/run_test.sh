

set -ex



pip check
pytest --cov=arrow --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=99 -k "not parse_tz_name_zzz"
exit 0
