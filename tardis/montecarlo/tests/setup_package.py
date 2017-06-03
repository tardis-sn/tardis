#setting the right include
from setuptools import Extension
import numpy as np
import os

from glob import glob


def get_extensions():
    sources = ['tardis/montecarlo/tests/packet_test_wrapper.pyx']
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), '..', 'src', '*.c'))]
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), '..', 'src',
                     'randomkit', '*.c'))]

    return [Extension('tardis.montecarlo.tests.packet_test_wrapper', sources,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    np.get_include()])]
