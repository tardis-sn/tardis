#setting the right include
from setuptools import Extension
import numpy as np
import os

randomkit_files = ['rk_isaac.c', 'rk_mt.c', 'rk_primitive.c','rk_sobol.c']

def get_extensions():
    return [Extension('tardis.montecarlo.montecarlo',
                      ['tardis/montecarlo/montecarlo.pyx',
                       'src/cmontecarlo.c'] +
                      [os.path.join('src/randomkit', fname)
                       for fname in randomkit_files],
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    np.get_include()])]
