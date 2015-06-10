import os
import random
import inspect
from ctypes import CDLL, c_double

import pytest
import numpy as np

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

tests = CDLL(test_path)

np.random.seed(1)

def generate_testcases_doubles(metafunc):
	from sys import float_info as floats
	MAX_FLOAT, MIN_FLOAT = floats[0], floats[3]
	for each in metafunc.fixturenames:
		if 'double_value' in inspect.getargspec(each):
			metafunc.parametrize("double_value", map(c_double, 
				np.random.uniform(MIN_FLOAT, MAX_FLOAT, size=10)))

def generate_testcases_integers(metafunc):
	MAX_INT, MIN_INT = 10000, -10000
	if 'int_value' in metafunc.fixturenames:
		metafunc.parametrize("int_value",
			np.random.randint(MIN_INT, MAX_INT, size=10))

def generate_testcases_unsigned_integers(metafunc):
	MAX_INT = 10000
	if 'unsigned_int_value' in metafunc.params:
		metafunc.parametrize("unsigned_int_value",
			np.random.randint(0, MAX_INT, size=10))


# Testing functions with float(C double) valued parameters
@pytest.fixture
def test_rpacket_get_nu(double_value):
	assert tests.test_rpacket_get_nu(request.params)

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