# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from astropy_helpers.setup_helpers import get_distutils_option
from setuptools import Extension


def get_package_data():
    return {
        "tardis.iip_plasma.tests": ["data/*.dat", "data/*.yml", "data/*.h5"]
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
    sources = ["tardis/iip_plasma/properties/util/macro_atom.pyx"]
    return [
        Extension(
            "tardis.iip_plasma.properties.util.macro_atom",
            sources,
            include_dirs=[np.get_include()],
            extra_compile_args=compile_args,
            extra_link_args=link_args,
        )
    ]
