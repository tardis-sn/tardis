

set -ex



pip check
python -c "import importlib; importlib.import_module('sphinx-jsonschema')"
exit 0
