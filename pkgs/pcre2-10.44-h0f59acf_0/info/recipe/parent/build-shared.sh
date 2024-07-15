#!/bin/bash
set -euxo pipefail

pushd build_cmake
ninja install

# Get static libs out of the way for now
mkdir -p ${SRC_DIR}/static_libs_for_cf
mv ${PREFIX}/lib/libpcre2*.a ${SRC_DIR}/static_libs_for_cf/
