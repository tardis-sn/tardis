#!/usr/bin/env bash

wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b

export PATH=/home/vagrant/miniconda/bin:$PATH
conda update --yes conda

conda create --yes -n tardis --file /vagrant/pip-requirements pip

echo "export PATH=/home/vagrant/miniconda/bin:$PATH" >> .bashrc

wget https://www.dropbox.com/s/svvyr5i7m8ouzdt/tardis_example.tar.gz
tar zxvf tardis_example.tar.gz