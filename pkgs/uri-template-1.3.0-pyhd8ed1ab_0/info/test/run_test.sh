

set -ex



pip check
pip list
pip list | grep -iE "uri-template\s+1\.3\.0"
mypy -p uri_template
cd src && coverage run --source=uri_template test.py && coverage report -m --fail-under=91
exit 0
