# Licensed under a 3-clause BSD style license - see LICENSE.rst
from setuptools import Extension
from astropy_helpers.distutils_helpers import get_distutils_option
import numpy as np


def get_package_data():
    return {
        "tardis.plasma.tests": [
            "data/*.dat",
            "data/*.txt",
            "data/*.yml",
            "data/*.h5",
            "data/*.dot",
            "data/*.tex",
        ]
    }


if (
    get_distutils_option("with_openmp", ["build", "install", "develop"])
    is not None
):
    compile_args = [
        "-fopenmp",
        "-W",
        "-Wall",
        "-Wmissing-prototypes",
        "-std=c99",
    ]
    link_args = ["-fopenmp"]
    define_macros = [("WITHOPENMP", None)]
else:
    compile_args = ["-W", "-Wall", "-Wmissing-prototypes", "-std=c99"]
    link_args = []
    define_macros = []


def get_extensions():
    sources = []
    return []
