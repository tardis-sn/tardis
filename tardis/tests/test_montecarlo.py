import os
from ctypes import *

from tardis import __path__ as path
from tardis import montecarlo

test_path = os.path.join(path[0], 'montecarlo', 'tests', 'test_cmontecarlo.so')
dl = CDLL(test_path)

def test():
	assert dl.testing_rpacket_get_nu() == True