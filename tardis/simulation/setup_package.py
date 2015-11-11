# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import Extension
from astropy_helpers.setup_helpers import get_distutils_option
import numpy as np

def get_package_data():
    return {'tardis.simulation.tests':['data/*.h5']}


