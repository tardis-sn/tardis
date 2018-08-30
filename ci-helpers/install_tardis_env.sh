#!/usr/bin/env bash
cd $TRAVIS_BUILD_DIR
if test -e $HOME/miniconda/envs/tardis; then
    echo "TARDIS env already installed.";
    # Also check for tardis_env27.yml change
else
   conda env create -f tardis_env27.yml
fi

source activate tardis