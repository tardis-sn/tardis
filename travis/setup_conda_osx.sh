#!/bin/bash

# Workaround for https://github.com/travis-ci/travis-ci/issues/6307, which
# caused the following error on MacOS X workers:
#
# /Users/travis/build.sh: line 109: shell_session_update: command not found
#
rvm get head

# Install conda
# http://conda.pydata.org/docs/travis.html#the-travis-yml-file
echo ls -la $HOME/miniconda
ls -la $HOME/miniconda
if [ ! "$(ls -A $HOME/miniconda)" ]; then
    echo "Cache is empty, installing miniconda and test environment."
    # Control will enter here if $HOME/miniconda is empty.
    rm -r $HOME/miniconda
    # Control will enter here if $HOME/miniconda doesn't exist.
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh

    export PATH="$HOME/miniconda/bin:$PATH"

    hash -r

    conda env create -f tardis_env27.yml
else
    export PATH="$HOME/miniconda/bin:$PATH"
fi

source activate tardis
