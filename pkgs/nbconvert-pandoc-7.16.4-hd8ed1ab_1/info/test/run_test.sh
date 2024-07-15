

set -ex



pip check
jupyter nbconvert --version
jupyter dejavu --version
jupyter nbconvert --help
jupyter dejavu --help
python check_pandoc.py ">=2.9.2" "<4.0.0"
exit 0
