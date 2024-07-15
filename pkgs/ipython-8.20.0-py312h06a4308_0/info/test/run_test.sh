

set -ex



pip check
pygmentize -L | grep ipython
ipython -h
ipython3 -h
exit 0
