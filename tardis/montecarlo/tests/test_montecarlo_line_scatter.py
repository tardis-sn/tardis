import os
from ctypes import CDLL, c_double, c_longdouble

import pytest

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

tests = CDLL(test_path)

def test_compute_distance2boundary():
	tests.test_compute_distance2boundary.restype = c_longdouble

	tests.init_rpacket()
	tests.init_storage_model()
	assert tests.test_compute_distance2boundary() == 0.5

# def test_montecarlo_line_scatter():
#	assert tests.test_montecarlo_line_scatter()