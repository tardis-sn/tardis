#!/usr/bin/env bash
cd $TRAVIS_BUILD_DIR
if test -e $HOME/miniconda/envs/tardis; then
    echo "TARDIS env already installed.";
    # Also check for tardis_env27.yml change
else
   conda env create -f tardis_env27.yml
   #trouble with building due to segfault at cython (https://github.com/cython/cython/issues/2199)
   #remove if we can get normal cython through conda
   source activate tardis
   conda uninstall -y cython
   git clone https://github.com/cython/cython
   cd cython
   git checkout c485b1b77264c3c75d090a3c526de24966830d42
   CFLAGS="$CFLAGS -D CYTHON_CLINE_IN_TRACEBACK=0" python setup.py install
   cd ..
fi

source activate tardis
