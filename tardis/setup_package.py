#setting the right include
from setuptools import Extension
import numpy as np

from astropy_helpers.setup_helpers import get_distutils_build_option

randomkit_files = ['tardis/randomkit/rk_isaac.c', 'tardis/randomkit/rk_mt.c',
                   'tardis/randomkit/rk_primitive.c',
                   'tardis/randomkit/rk_sobol.c']

def get_extensions():
    return [Extension('tardis.montecarlo',
                      ['tardis/montecarlo.pyx', 'tardis/cmontecarlo.c'] +
                      randomkit_files,
                      include_dirs=['tardis/randomkit', np.get_include()],
                      extra_compile_args=['-fopenmp'],
                      extra_link_args=['-lgomp']
                      )]

def _get_compile_link_args():
    no_openmp = get_distutils_build_option('no_openmp')

    if no_openmp is None or no_openmp==False:
        extra_compile_args = ['-fopenmp']
        extra_link_args = ['-lgomp']

    else:
        extra_compile_args = []
        extra_link_args = []

    return extra_compile_args, extra_link_args
