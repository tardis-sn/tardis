#!/bin/bash

cat <<EOF > setup.py
import os
from setuptools import setup


version = os.environ['PKG_VERSION']

setup(name="dataclasses",
      version=version,
      platforms=["any"],
      python_requires='>=3.7',
      description='Empty dataclasses shim to make it visible for pkg_resources under Python 3.7+'
)
EOF

$PYTHON -m pip install --no-deps --ignore-installed .