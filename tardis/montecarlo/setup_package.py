#setting the right include
from setuptools import Extension
import numpy as np
import os
from astropy_helpers.setup_helpers import get_distutils_option
from Cython.Build import cythonize

from glob import glob

if get_distutils_option('with_openmp', ['build', 'install', 'develop']) is not None:
    compile_args = ['-fopenmp', '-W', '-Wall', '-Wmissing-prototypes', '-std=c99']
    link_args = ['-fopenmp']
    define_macros = [('WITHOPENMP', None)]
else:
    compile_args = ['-W', '-Wall', '-Wmissing-prototypes', '-std=c99']
    link_args = []
    define_macros = []

if get_distutils_option('with_vpacket_logging', ['build', 'install', 'develop']) is not None:
    define_macros.append(('WITH_VPACKET_LOGGING', None))

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
    #tests = [os.path.relpath(fname) for fname in os.wakk(
    #    os.path.join(os.path.dirname(__file__), 'src', '*.c'))
    #        if os.path.basename(fname).startswith('test')]
    test_sources = [
            os.path.relpath(os.path.join(r,fname)) for r,p,fnames in os.walk(os.path.dirname(__file__))
            for fname in fnames
            if os.path.basename(fname).endswith('.c')
            ]
    test_deps = [
            os.path.relpath(os.path.join(r,fname)) for r,p,fnames in os.walk(os.path.dirname(__file__))
            for fname in fnames
            if os.path.basename(fname).endswith('.h')
            ]

    return cythonize(
            Extension('tardis.montecarlo.montecarlo', sources,
                include_dirs=['tardis/montecarlo/src',
                    'tardis/montecarlo/src/randomkit',
                    'numpy'],
                depends=deps,
                extra_compile_args=compile_args,
                extra_link_args=link_args,
                define_macros=define_macros)
            ) + [
# Test extension
            Extension('tardis.montecarlo.test_montecarlo', test_sources,
                include_dirs=['tardis/montecarlo/src',
                    'tardis/montecarlo/src/randomkit',
                    'numpy'],
                depends=test_deps,
                extra_compile_args=compile_args,
                extra_link_args=link_args,
#                library_dirs=['tardis/montecarlo'],
#                libraries= ['montecarlo'],
                define_macros=define_macros),
            ]

def get_package_data():
    return {'tardis.montecarlo.tests':['data/*.npy']}
