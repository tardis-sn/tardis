import os
from sys import float_info as floats
import random
from ctypes import CDLL, c_double

import pytest
import numpy as np

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')
MAX_FLOAT, MIN_FLOAT = floats[0], floats[3]
MAX_INT = 10000
MIN_INT = -10000

tests = CDLL(test_path)
SIZE = 10


# Testing functions with float(C double) valued parameters
np.random.seed(1)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_nu(double_value):
	assert tests.test_rpacket_get_nu(double_value)

np.random.seed(2)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_mu(double_value):
	assert tests.test_rpacket_get_mu(double_value)

np.random.seed(3)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_energy(double_value):
	assert tests.test_rpacket_get_energy(double_value)

np.random.seed(4)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_r(double_value):
	assert tests.test_rpacket_get_r(double_value)

np.random.seed(5)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_tau_event(double_value):
	assert tests.test_rpacket_get_tau_event(double_value)

np.random.seed(6)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_nu_line(double_value):
	assert tests.test_rpacket_get_nu_line(double_value)

np.random.seed(7)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_nu_line(double_value):
	assert tests.test_rpacket_get_d_boundary(double_value)

np.random.seed(8)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_nu_line(double_value):
	assert tests.test_rpacket_get_d_electron(double_value)

np.random.seed(9)
@pytest.mark.parametrize(
	"double_value", map(c_double, np.random.uniform(
						MIN_FLOAT, MAX_FLOAT, size=SIZE))
	)
def test_rpacket_get_nu_line(double_value):
	assert tests.test_rpacket_get_d_line(double_value)


# Testing functions with Integer valued parameters
np.random.seed(10)
@pytest.mark.parametrize(
	"unsigned_int_value", np.random.randint(0, MAX_INT, size=SIZE)
	)
def test_rpacket_get_current_shell_id(unsigned_int_value):
	assert tests.test_rpacket_get_current_shell_id(unsigned_int_value)

np.random.seed(11)
@pytest.mark.parametrize(
	"nt_value", np.random.randint(MIN_INT, MAX_INT, size=SIZE)
	)
def test_rpacket_get_next_line_id(int_value):
	assert tests.test_rpacket_get_next_line_id(int_value)

np.random.seed(12)
@pytest.mark.parametrize(
	"int_value", np.random.randint(MIN_INT, MAX_INT, size=SIZE)
	)
def test_rpacket_get_next_line_id(int_value):
	assert tests.test_rpacket_get_recently_crossed_boundary(int_value)

np.random.seed(13)
@pytest.mark.parametrize(
	"int_value", np.random.randint(MIN_INT, MAX_INT, size=SIZE)
	)
def test_rpacket_get_next_line_id(int_value):
	assert tests.test_rpacket_get_virtual_packet_flag(int_value)

np.random.seed(14)
@pytest.mark.parametrize(
	"int_value", np.random.randint(MIN_INT, MAX_INT, size=SIZE)
	)
def test_rpacket_get_next_line_id(int_value):
	assert tests.test_rpacket_get_virtual_packet(int_value)

np.random.seed(15)
@pytest.mark.parametrize(
	"int_value", np.random.randint(MIN_INT, MAX_INT, size=SIZE)
	)
def test_rpacket_get_next_line_id(int_value):
	assert tests.test_rpacket_get_next_line_id(int_value)


# Testing functions without any parameters
def test_rpacket_get_last_line():
	assert tests.test_rpacket_get_last_line()

def test_rpacket_get_close_line():
	assert tests.test_rpacket_get_close_line()

def test_rpacket_get_status():
	assert tests.test_rpacket_get_status()