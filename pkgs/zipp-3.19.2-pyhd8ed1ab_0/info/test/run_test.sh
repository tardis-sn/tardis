

set -ex



python -m unittest tests/test_path.py
pip check
exit 0
