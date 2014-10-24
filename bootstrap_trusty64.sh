#!/usr/bin/env bash

apt-get update
apt-get install -y python-virtualenv python-numpy python-pandas python-scipy python-h5py python-yaml python-astropy cython
chown -R
virtualenv --system-site-packages tardis.ve
source tardis.ve/bin/activate
pip install -r /vagrant/pip-requirements
