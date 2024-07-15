#!/bin/bash

source ${RECIPE_DIR}/setup_compiler.sh

if [[ "${PKG_NAME}" == "gcc" ]]; then
  for tool in cc cpp gcc gcc-ar gcc-nm gcc-ranlib gcov gcov-dump gcov-tool; do
    ln -sf ${PREFIX}/bin/${triplet}-${tool} ${PREFIX}/bin/${tool}
  done
elif [[ "${PKG_NAME}" == "gxx" ]]; then
  ln -sf ${PREFIX}/bin/${triplet}-g++ ${PREFIX}/bin/g++
  ln -sf ${PREFIX}/bin/${triplet}-c++ ${PREFIX}/bin/c++
elif [[ "${PKG_NAME}" == "gfortran" ]]; then
  ln -sf ${PREFIX}/bin/${triplet}-gfortran ${PREFIX}/bin/gfortran
fi
