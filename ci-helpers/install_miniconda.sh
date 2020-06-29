#!/usr/bin/env bash
export PYTHONIOENCODING=UTF8
if test -e $HOME/miniconda/bin; then
    echo "miniconda already installed.";
    export PATH=$HOME/miniconda/bin:$PATH
    hash -r
    #conda update --yes conda

    # install mamba
    conda install mamba -c conda-forge

else
    wget $MINICONDA_URL -O miniconda.sh
    chmod +x miniconda.sh
    rm -r $HOME/miniconda
    bash miniconda.sh -b -p $HOME/miniconda

    export PATH=$HOME/miniconda/bin:$PATH
    hash -r
    conda update --yes conda

    # install mamba
    conda install mamba -c conda-forge

fi

source $HOME/miniconda/etc/profile.d/conda.sh