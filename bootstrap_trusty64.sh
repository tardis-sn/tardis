#!/usr/bin/env bash

apt-get update
apt-get install -y python-virtualenv python-numpy python-pandas python-scipy python-h5py python-yaml ipython python-matplotlib cython git
virtualenv --system-site-packages tardis.ve
chmod -R a+rwX tardis.ve
source tardis.ve/bin/activate
pip install -r /vagrant/pip-requirements
