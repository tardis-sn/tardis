#!/usr/bin/env python

import sys

from setuptools_scm import version_from_scm

version = version_from_scm(".").tag.public

print(version)
sys.exit(0)
