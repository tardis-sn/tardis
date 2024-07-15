

set -ex



pip check
conda install -y numpy
python test-imports-requiring-numpy.py
conda install -y matplotlib-base
python test-imports-requiring-matplotlib.py
conda install -y pandas
python test-imports-requiring-pandas.py
exit 0
