#!/bin/bash

cat <<EOF > conda_mpl_config.cfg
[directories]
basedirlist = $PREFIX

[packages]
tests = False
toolkit_tests = False
sample_data = False

[libs]
system_freetype = True

EOF

cat conda_mpl_config.cfg
sed -i.bak "s|/usr/local|${PREFIX}|" setupext.py

export MPLSETUPCFG=conda_mpl_config.cfg

$PYTHON -m pip install --no-deps --no-build-isolation -vv .
