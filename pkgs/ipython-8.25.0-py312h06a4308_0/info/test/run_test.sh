

set -ex



pip check
pygmentize -L | grep ipython
ipython -h
ipython3 -h
export MIGRATING=false
exit 0
