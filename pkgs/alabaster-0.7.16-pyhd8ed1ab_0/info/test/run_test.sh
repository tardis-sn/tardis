

set -ex



pip check
python -c "import alabaster; print(alabaster.get_path())"
exit 0
