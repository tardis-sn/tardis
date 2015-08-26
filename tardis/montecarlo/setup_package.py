#setting the right include
from setuptools import Extension
import numpy as np
import os
from astropy_helpers.setup_helpers import get_distutils_option

from glob import glob

if get_distutils_option('with_openmp', ['build', 'install', 'develop']) is not None:
    compile_args = ['-fopenmp', '-W', '-Wall', '-Wmissing-prototypes', '-std=c99']
    link_args = ['-fopenmp']
    define_macros = [('WITHOPENMP', None)]
else:
    compile_args = ['-W', '-Wall', '-Wmissing-prototypes', '-std=c99']
    link_args = []
    define_macros = []

def get_extensions():
    sources = ['tardis/montecarlo/montecarlo.pyx']
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src', '*.c'))
                if not os.path.basename(fname).startswith('test')]
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.c'))]
    deps = [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src', '*.h'))]
    deps += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.h'))]

    return [Extension('tardis.montecarlo.montecarlo', sources,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    np.get_include()],
                      depends=deps,
                      extra_compile_args=compile_args,
                      extra_link_args=link_args,
                      define_macros=define_macros)]
