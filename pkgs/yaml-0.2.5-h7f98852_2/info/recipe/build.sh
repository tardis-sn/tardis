#!/bin/bash

./bootstrap
./configure --prefix="${PREFIX}"
make
if [[ "${CONDA_BUILD_CROSS_COMPILATION}" != "1" ]]; then
make check
fi
make install
