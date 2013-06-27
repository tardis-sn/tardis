__author__ = 'michi'

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import glob
import os

import sys
#from setuptools import setup, find_packages
import numpy as np

# Set affiliated package-specific settings
PACKAGENAME = 'tardis'
DESCRIPTION = 'TARDIS'
LONG_DESCRIPTION = ''
AUTHOR = ''
AUTHOR_EMAIL = ''
LICENSE = 'BSD'
URL = 'http://astropy.org'

#version should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
version = '0.0.dev'

# Indicates if this version is a release version
release = 'dev' not in version

#Define the extension randomkit
cflags = ['-fPIC']
cfiles = ['tardis/randomkit/rk_isaac.c', 'tardis/randomkit/rk_mt.c', 'tardis/randomkit/rk_primitive.c',
          'tardis/randomkit/rk_sobol.c']

randmomkit = Extension('randomkit',
    cfiles,
    extra_compile_args=cflags)

montecarlo_multizone = Extension('tardis.montecarlo_multizone',
    ['tardis/montecarlo_multizone.pyx'] + cfiles)

macro_atom = Extension('tardis.macro_atom',
    ['tardis/macro_atom.pyx'])

#test_cython = Extension('tardis.test_cython',
#    ['tardis/test_cython.pyx'])

# Define all packages  and modules
packages = ['tardis', 'tardis.tests']
package_dir = {'tardis': 'tardis'}
package_data = {'tardis': ['data/*', 'tests/data/*']}

# Treat everything in scripts except README.rst as a script to be installed
scripts = glob.glob(os.path.join('scripts', '*'))
scripts.remove(os.path.join('scripts', 'README.rst'))

#Define the extension modules
extensions = [montecarlo_multizone, macro_atom]

setup(name=PACKAGENAME,
    version=version,
    description=DESCRIPTION,
    packages=packages,
    package_dir=package_dir,
    package_data=package_data,
    ext_modules=extensions,
    include_dirs=[np.get_include()],
    scripts=scripts,
    requires=['astropy', 'numpy', 'scipy'],
    install_requires=['astropy'],
    provides=[PACKAGENAME],
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    url=URL,
    long_description=LONG_DESCRIPTION,
    cmdclass={'build_ext': build_ext}
    #zip_safe=False,
    #use_2to3=True
)



