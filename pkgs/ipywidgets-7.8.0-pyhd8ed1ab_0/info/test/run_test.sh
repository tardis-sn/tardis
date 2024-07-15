

set -ex



pip check
pytest -vv --pyargs ipywidgets --cov=ipywidgets --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=91
exit 0
