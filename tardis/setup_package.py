#setting the right include
from setuptools import Extension
import numpy as np

randomkit_files = ['tardis/randomkit/rk_isaac.c', 'tardis/randomkit/rk_mt.c',
                   'tardis/randomkit/rk_primitive.c',
                   'tardis/randomkit/rk_sobol.c']

def get_extensions():
    return [Extension('tardis.montecarlo',
                        ['tardis/montecarlo.pyx', 'tardis/cmontecarlo.c'] +
                        randomkit_files,
                        include_dirs=['tardis/randomkit', np.get_include()]),
            Extension('tardis.tests.montecarlo_test_wrappers',
                      ['tardis/tests/montecarlo_test_wrappers.pyx',
                       'tardis/cmontecarlo.c'] + randomkit_files,
                      include_dirs=['tardis/randomkit', np.get_include()])]

