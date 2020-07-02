#setting the right include
from tardis import __path__ as TARDIS_PATH
from setuptools import Extension
import os
from astropy_helpers.distutils_helpers import get_distutils_option
# from Cython.Build import cythonize

import yaml

from glob import glob

TARDIS_PATH = TARDIS_PATH[0]

if get_distutils_option('with_openmp', ['build', 'install', 'develop']) is not None:
    compile_args = ['-fopenmp', '-W', '-Wall', '-Wmissing-prototypes', '-std=c99']
    link_args = ['-fopenmp']
    define_macros = [('WITHOPENMP', None)]
else:
    compile_args = ['-W', '-Wall', '-Wmissing-prototypes', '-std=c99']
    link_args = []
    define_macros = []


vpacket_config_file_path = os.path.join(TARDIS_PATH, 'data', 'vpacket_config.yml')
if get_distutils_option('with_vpacket_logging', ['build', 'install', 'develop']) is not None:
    define_macros.append(('WITH_VPACKET_LOGGING', None))
    vpacket_config = {'vpacket_logging': True}
else:
    vpacket_config = {'vpacket_logging': False}

yaml.dump(vpacket_config, open(vpacket_config_file_path, "w"), default_flow_style=False)

# def get_extensions():
#     sources = ['tardis/montecarlo/montecarlo.pyx']
#     sources += [os.path.relpath(fname) for fname in glob(
#         os.path.join(os.path.dirname(__file__), 'src', '*.c'))]
#     sources += [os.path.relpath(fname) for fname in glob(
#         os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.c'))]
#     deps = [os.path.relpath(fname) for fname in glob(
#         os.path.join(os.path.dirname(__file__), 'src', '*.h'))]
#     deps += [os.path.relpath(fname) for fname in glob(
#         os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.h'))]

    # return cythonize(
    #         Extension('tardis.montecarlo.montecarlo', sources,
    #             include_dirs=['tardis/montecarlo/src',
    #                 'tardis/montecarlo/src/randomkit',
    #                 'numpy'],
    #             depends=deps,
    #             extra_compile_args=compile_args,
    #             extra_link_args=link_args,
    #             define_macros=define_macros)
    #         )


def get_package_data():
    return {'tardis.montecarlo.tests':['data/*.npy', 'data/*.hdf']}
