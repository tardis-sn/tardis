#!/bin/bash
cat > site.cfg <<EOF
[mkl]
library_dirs = $PREFIX/lib
include_dirs = $PREFIX/include
libraries = mkl_rt
EOF
echo "#####################"
echo "site.cfg file written"

$PYTHON -m pip install . --no-deps -vv
