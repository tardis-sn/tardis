#setting the right include
from setuptools import Extension
import os
from glob import glob

import numpy as np

def get_extensions():
    sources = ['tardis/montecarlo/montecarlo.pyx']
    test_sources = []
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src', '*.c'))]
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.c'))]
    test_sources += [os.path.relpath(fname) for fname in \
        glob(os.path.join(os.path.dirname(__file__), 'tests', '*.c'))]

    return [Extension('tardis.montecarlo.montecarlo', sources,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    np.get_include()]),
            Extension('tardis.montecarlo.tests.test_cmontecarlo', test_sources,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    np.get_include()])
            ]