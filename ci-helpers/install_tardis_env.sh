#!/usr/bin/env bash
cd $TRAVIS_BUILD_DIR
if test -e  $HOME/miniconda/envs/tardis; then
    echo "TARDIS env already installed.";
    # Also check for tardis_env3.yml change
else
   conda env create -f tardis_env3.yml
   #trouble with building due to segfault at cython (https://github.com/cython/cython/issues/2199)
   #remove if we can get normal cython through conda
   conda activate tardis
fi

conda activate tardis
