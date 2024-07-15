#!/bin/bash

export BLOSC_DIR=$PREFIX
export BLOSC2_DIR=$PREFIX
export BZIP2_DIR=$PREFIX
export HDF5_DIR=$PREFIX
export LZO_DIR=$PREFIX
export COPY_DLLS=FALSE
export PYTABLES_NO_BLOSC2_WHEEL=TRUE

# Remove the pre-cythonized files which may not be compatible.
rm -f tables/*.c

$PYTHON -m pip install --no-deps --no-cache-dir --ignore-installed .
