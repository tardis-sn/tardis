

set -ex



pip check
mypy -m async_lru
pytest -vv --asyncio-mode=auto --cov=async_lru --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=89
exit 0
