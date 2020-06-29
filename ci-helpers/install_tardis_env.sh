#!/usr/bin/env bash
cd $TARDIS_BUILD_DIR

# install mamba
if test -e $HOME/miniconda/bin/mamba; then
    echo "Mamba installed Succesfuly"

else
    echo "##vso[task.loginfo type=debug;] Mamba not installed"
    conda install mamba -c conda-forge -y
    echo "##vso[task.loginfo type=debug;] Installed Mamba to correct location"
fi

if test -e  $HOME/miniconda/envs/tardis; then
    echo "TARDIS env already installed.";
    # Also check for tardis_env3.yml change
else
   # conda env create -f tardis_env3.yml
   mamba env create -f tardis_env3.yml
   #trouble with building due to segfault at cython (https://github.com/cython/cython/issues/2199)
   #remove if we can get normal cython through conda
fi
