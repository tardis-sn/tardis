

set -ex



pip check
pip list | grep -iE nbformat
pip list | grep -iE 'nbformat\s*5\.10\.4'
python -c "v = __import__('nbformat').__version__; print(v); assert v == '5.10.4'"
jupyter trust --version
jupyter-trust --help
pytest -vv --cov=nbformat --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=76
exit 0
