#!/usr/bin/env python

import sys
from datetime import date
from setuptools_scm import version_from_scm
from setuptools_scm.version import guess_next_date_ver

scm_version = version_from_scm(".").tag.public
scm_version = guess_next_date_ver(scm_version)
scm_version = [ int(i) for i in scm_version.split(".") ]

build = scm_version[3]
release = scm_version[:3]
iso_date = date(*release).isoformat()
iso_date = iso_date.replace("-",".")

version = f"{iso_date}.{str(build)}" if build > 0 else iso_date
print(version)

sys.exit(0)
