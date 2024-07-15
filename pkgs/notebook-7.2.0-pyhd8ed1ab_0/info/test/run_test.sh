

set -ex



pip check
jupyter notebook -h
jupyter-notebook -h
jupyter labextension list
jupyter labextension list 1>labextensions 2>&1
grep -iE "@jupyter-notebook/lab-extension.*OK" labextensions
exit 0
