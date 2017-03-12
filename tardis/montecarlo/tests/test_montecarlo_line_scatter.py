import os
from ctypes import CDLL, c_double

import pytest

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

tests = CDLL(test_path)

tests.init_rpacket()
tests.init_storage_model()

def test_compute_distance2boundary():
	tests.test_compute_distance2boundary.restype = c_double
	assert tests.test_compute_distance2boundary() == 7572770063489.144

# def test_montecarlo_line_scatter():
#	assert tests.test_montecarlo_line_scatter()