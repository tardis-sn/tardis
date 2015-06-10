import os
import random
import inspect
from ctypes import CDLL, c_double

import pytest
import numpy as np

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

tests = CDLL(test_path)

np.random.seed(26061963)


@pytest.fixture(params=np.random.random(size=10))
def random_double_value(request):
	return request.params



# Testing functions with float(C double) valued parameters

def test_rpacket_get_nu(random_double_value):
	assert tests.test_rpacket_get_nu(random_double_value)

@pytest.fixture
def test_rpacket_get_mu(double_value):
	assert tests.test_rpacket_get_mu(double_value)

@pytest.fixture
def test_rpacket_get_energy(double_value):
	assert tests.test_rpacket_get_energy(double_value)

@pytest.fixture
def test_rpacket_get_r(double_value):
	assert tests.test_rpacket_get_r(double_value)

@pytest.fixture
def test_rpacket_get_tau_event(double_value):
	assert tests.test_rpacket_get_tau_event(double_value)

@pytest.fixture
def test_rpacket_get_nu_line(double_value):
	assert tests.test_rpacket_get_nu_line(double_value)

@pytest.fixture
def test_rpacket_get_d_boundary(double_value):
	assert tests.test_rpacket_get_d_boundary(double_value)

@pytest.fixture
def test_rpacket_get_d_electron(double_value):
	assert tests.test_rpacket_get_d_electron(double_value)

@pytest.fixture
def test_rpacket_get_d_line(double_value):
	assert tests.test_rpacket_get_d_line(double_value)


# Testing functions with Unsigned Integer valued parameters
@pytest.fixture
def test_rpacket_get_current_shell_id(unsigned_int_value):
	assert tests.test_rpacket_get_current_shell_id(unsigned_int_value)

@pytest.fixture
def test_rpacket_get_next_line_id(unsigned_int_value):
	assert tests.test_rpacket_get_next_line_id(unsigned_int_value)


# Testing functions with Integer valued parameters
@pytest.fixture
def test_rpacket_get_recently_crossed_boundary(int_value):
	assert tests.test_rpacket_get_recently_crossed_boundary(int_value)

@pytest.fixture
def test_rpacket_get_virtual_packet_flag(int_value):
	assert tests.test_rpacket_get_virtual_packet_flag(int_value)

@pytest.fixture
def test_rpacket_get_virtual_packet(int_value):
	assert tests.test_rpacket_get_virtual_packet(int_value)

@pytest.fixture
def test_rpacket_get_next_shell_id(int_value):
	assert tests.test_rpacket_get_next_shell_id(int_value)


# Testing functions without any parameters
def test_rpacket_get_last_line():
	assert tests.test_rpacket_get_last_line()

def test_rpacket_get_close_line():
	assert tests.test_rpacket_get_close_line()

def test_rpacket_get_status():
	assert tests.test_rpacket_get_status()