# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import glob
import os
import sys

import setuptools_bootstrap
from setuptools import setup, Extension

import numpy as np
#A dirty hack to get around some early import/configurations ambiguities
if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
builtins._ASTROPY_SETUP_ = True

import astropy
from astropy.setup_helpers import (register_commands, adjust_compiler,
                                   get_debug_option)
from astropy.version_helpers import get_git_devstr, generate_version_py


PACKAGENAME = 'tardis'
DESCRIPTION = 'TARDIS - Temperature And Radiative Diffusion In Supernovae'
LONG_DESCRIPTION = ''
AUTHOR = 'TARDIS Team'
AUTHOR_EMAIL = 'wkerzendorf@gmail.com'
LICENSE = 'BSD (3 clause)'
URL = 'http://tardis.readthedocs.org/'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.9'

# Indicates if this version is a release version
RELEASE = 'dev' not in VERSION

if not RELEASE:
    VERSION += get_git_devstr(False)

# Populate the dict of setup command overrides; this should be done before
# invoking any other functionality from distutils since it can potentially
# modify distutils' behavior.
cmdclassd = register_commands(PACKAGENAME, VERSION, RELEASE)

# Adjust the compiler in case the default on this platform is to use a
# broken one.
adjust_compiler(PACKAGENAME)

# Freeze build information in version.py
generate_version_py(PACKAGENAME, VERSION, RELEASE, get_debug_option())

# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
           if fname != 'README.rst']


# Use the find_packages tool to locate all packages and modules

try:

    from astropy.setup_helpers import get_package_info

    # Get configuration information from all of the various subpackages.
    # See the docstring for setup_helpers.update_package_files for more
    # details.
    package_info = get_package_info(PACKAGENAME)

    # Add the project-global data
    package_info['package_data'][PACKAGENAME] = ['data/*']

except ImportError: # compatibility with Astropy 0.2 - can be removed in cases
                    # where Astropy 0.2 is no longer supported

    from setuptools import find_packages
    from astropy.setup_helpers import filter_packages, update_package_files

    package_info = {}

    # Use the find_packages tool to locate all packages and modules
    package_info['packages'] = filter_packages(find_packages())

    # Additional C extensions that are not Cython-based should be added here.
    package_info['ext_modules'] = []

    # A dictionary to keep track of all package data to install
    package_info['package_data'] = {PACKAGENAME: ['data/*']}

    # A dictionary to keep track of extra packagedir mappings
    package_info['package_dir'] = {}

    # Update extensions, package_data, packagenames and package_dirs from
    # any sub-packages that define their own extension modules and package
    # data.  See the docstring for setup_helpers.update_package_files for
    # more details.
    update_package_files(PACKAGENAME, package_info['ext_modules'],
                         package_info['package_data'], package_info['packages'],
                         package_info['package_dir'])


# Additional C extensions that are not Cython-based should be added here.


#extensions = [Extension('tardis.montecarlo',
#                        ['tardis/montecarlo.pyx'] + randomkit_files)]



setup(name=PACKAGENAME+'-sn',
      version=VERSION,
      description=DESCRIPTION,
      scripts=scripts,
      include_dirs=[np.get_include()],
      install_requires=['astropy>=0.3', 'numpy', 'scipy', 'h5py', 'pandas>=0.12', 'pyyaml'],
      #setup_requires=['astropy>=0.2.4', 'numpy'],
      provides=[PACKAGENAME],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url=URL,
      long_description=LONG_DESCRIPTION,
      cmdclass=cmdclassd,
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: C',
          'Programming Language :: Cython',
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: Implementation :: CPython',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics'
      ],
      zip_safe=False,
      use_2to3=True,
      **package_info
)
