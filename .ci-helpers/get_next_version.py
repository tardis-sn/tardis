#!/usr/bin/env python

import sys
from setuptools_scm import version_from_scm
from setuptools_scm.version import guess_next_date_ver

version = version_from_scm(".").tag.public
version = guess_next_date_ver(version)

version = version.rstrip(".0") if version.endswith(".0") else version
version = version.split(".")

version[1:3] = ["0" + i if len(i) == 1 else i for i in version[1:3]]
version = ".".join(version)

print(version)

sys.exit(0)
