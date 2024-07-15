#!/bin/bash
set -euxo pipefail

mv ${SRC_DIR}/static_libs_for_cf/libpcre2*.a ${PREFIX}/lib/
