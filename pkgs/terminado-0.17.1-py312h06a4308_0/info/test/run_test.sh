

set -ex



pip check
pip list
pip list | grep -iE "terminado\s*0\.17\.1"
exit 0
