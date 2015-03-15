#setting the right include
from setuptools import Extension
import numpy as np
import os

from glob import glob


def get_extensions():
    sources = ['tardis/montecarlo/src/cmontecarlo.pyx']
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), '*.c'))]
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'randomkit', '*.c'))]

    return [Extension('tardis.montecarlo.src.montecarlo', sources,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    np.get_include()])]
