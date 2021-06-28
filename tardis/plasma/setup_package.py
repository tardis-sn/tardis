# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from setuptools import Extension

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
