#setting the right include
from setuptools import Extension
import numpy as np
import os

from glob import glob


def get_extensions():
    sources = ['tardis/montecarlo/montecarlo.pyx']
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src', '*.c'))]
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.c'))]

    sources_test = ['tardis/montecarlo/tests/test_cmontecarlo.c']
    sources_test += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src', '*.c'))]
    sources_test += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.c'))]

    return [Extension('tardis.montecarlo.montecarlo', sources,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    'tardis/montecarlo/tests',
                                    np.get_include()]),
	    Extension('tardis.montecarlo.tests.test_cmontecarlo', sources_test,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    'tardis/montecarlo/tests',
                                    np.get_include()])]
