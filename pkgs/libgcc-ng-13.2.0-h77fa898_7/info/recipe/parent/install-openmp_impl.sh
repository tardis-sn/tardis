#!/bin/bash

source ${RECIPE_DIR}/setup_compiler.sh
set -e -x

export CHOST="${triplet}"

mkdir -p ${PREFIX}/lib

pushd ${PREFIX}/lib/
ln -s libgomp.so.${libgomp_ver} libgomp.so.${libgomp_ver:0:1}
popd
