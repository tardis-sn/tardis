import os
from ctypes import CDLL, c_double

import pytest

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

tests = CDLL(test_path)

def test_montecarlo_line_scatter():
	assert tests.test_montecarlo_line_scatter()