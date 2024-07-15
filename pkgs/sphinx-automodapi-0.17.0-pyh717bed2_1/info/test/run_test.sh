

set -ex



python -m pip check sphinx-automodapi
python -m pip show sphinx-automodapi
python -m pytest -ra --pyargs sphinx_automodapi -k 'not cython'
exit 0
