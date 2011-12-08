#setup script for TARDIS software

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import glob

#version number
version = '0.01dev'

# Treat everything in scripts except README.rst as a script to be installed
scripts = glob.glob('scripts/*')
try:
	scripts.remove('scripts/README.rst')
except ValueError:
	pass


setup(name='tardis',
	  description='TARDIS Software - Time And Relative Diffusion in Supernovae',
	  author='The TARDIS collaboration',
      version=version,
      packages=['tardis'],
      package_data={'tardis': ['data/*']},
      scripts=scripts
      )
      
