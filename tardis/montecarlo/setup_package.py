#setting the right include
from setuptools import Extension
import numpy as np
import os
from astropy_helpers.setup_helpers import get_distutils_option

from glob import glob

if get_distutils_option('with_openmp', ['build', 'install']) is not None:
    compile_args = ['-fopenmp']
    link_args = ['-fopenmp']
    define_macros = [('WITHOPENMP', None)]
else:
    compile_args = []
    link_args = []
    define_macros = []

def get_extensions():
    sources = ['tardis/montecarlo/montecarlo.pyx']
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src', '*.c'))]
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.c'))]

    return [Extension('tardis.montecarlo.montecarlo', sources,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    np.get_include()],
                      extra_compile_args=compile_args,
                      extra_link_args=link_args,
                      define_macros=define_macros)]
