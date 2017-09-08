#!/bin/bash

# Install conda
# http://conda.pydata.org/docs/travis.html#the-travis-yml-file
echo ls -la $HOME/miniconda
ls -la $HOME/miniconda
if [ ! "$(ls -A $HOME/miniconda)" ]; then
    echo "Cache is empty, installing miniconda and test environment."
    # Control will enter here if $HOME/miniconda is empty.
    rm -r $HOME/miniconda
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda

    export PATH="$HOME/miniconda/bin:$PATH"

    hash -r

    conda env create -f tardis_env27.yml
else
    export PATH="$HOME/miniconda/bin:$PATH"
fi

source activate tardis
