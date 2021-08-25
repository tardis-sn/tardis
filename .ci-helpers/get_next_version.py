#!/usr/bin/env python

import sys
from setuptools_scm import version_from_scm
from setuptools_scm.version import guess_next_date_ver

version = version_from_scm(".").tag.public
version = guess_next_date_ver(version)
version = version.rstrip(".0") if version.endswith(".0") else version

print(version)
sys.exit(0)
