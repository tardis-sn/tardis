

set -ex



python -c "from dataclasses import dataclass"
python -c "import pkg_resources as p; p.get_distribution('dataclasses')"
exit 0
