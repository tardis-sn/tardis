# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import Extension
import numpy as np

def get_package_data():
    return {'tardis.plasma.tests':['data/*.dat', 'data/*.yml']}

def get_extensions():
    sources = ['tardis/plasma/properties/util/macro_atom.pyx']
    return [Extension('tardis.plasma.properties.util.macro_atom', sources,
                      include_dirs=[np.get_include()])]