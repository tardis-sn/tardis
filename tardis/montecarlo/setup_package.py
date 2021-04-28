# setting the right include
from setuptools import Extension
from astropy_helpers.distutils_helpers import get_distutils_option
import numpy as np
from glob import glob


def get_package_data():
    return {"tardis.montecarlo.tests": ["data/*.npy", "data/*.hdf"]}

