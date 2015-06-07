import os
from sys import float_info as floats
import random
from ctypes import CDLL, c_double

import pytest
import numpy as np

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

tests = CDLL(test_path)
SIZE = 10

np.random.seed(1)
doubles = map(c_double, np.random.uniform(
			floats[3], floats[0], size=SIZE))
unsigned_integers = np.random.randint(0, 10000, size=SIZE)
integers = np.random.randint(-10000, 10000, size=SIZE)

# Testing functions with float(C double) valued parameters
@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_nu(double_value):
	assert tests.test_rpacket_get_nu(double_value)

@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_mu(double_value):
	assert tests.test_rpacket_get_mu(double_value)

@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_energy(double_value):
	assert tests.test_rpacket_get_energy(double_value)

@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_r(double_value):
	assert tests.test_rpacket_get_r(double_value)

@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_tau_event(double_value):
	assert tests.test_rpacket_get_tau_event(double_value)

@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_nu_line(double_value):
	assert tests.test_rpacket_get_nu_line(double_value)

@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_d_boundary(double_value):
	assert tests.test_rpacket_get_d_boundary(double_value)

@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_d_electron(double_value):
	assert tests.test_rpacket_get_d_electron(double_value)

@pytest.mark.parametrize(
	"double_value", doubles
	)
def test_rpacket_get_d_line(double_value):
	assert tests.test_rpacket_get_d_line(double_value)


# Testing functions with Unsigned Integer valued parameters
@pytest.mark.parametrize(
	"unsigned_int_value", unsigned_integers
	)
def test_rpacket_get_current_shell_id(unsigned_int_value):
	assert tests.test_rpacket_get_current_shell_id(unsigned_int_value)

@pytest.mark.parametrize(
	"unsigned_int_value", unsigned_integers
	)
def test_rpacket_get_next_line_id(unsigned_int_value):
	assert tests.test_rpacket_get_next_line_id(unsigned_int_value)


# Testing functions with Integer valued parameters
@pytest.mark.parametrize(
	"int_value", integers
	)
def test_rpacket_get_recently_crossed_boundary(int_value):
	assert tests.test_rpacket_get_recently_crossed_boundary(int_value)

@pytest.mark.parametrize(
	"int_value", integers
	)
def test_rpacket_get_virtual_packet_flag(int_value):
	assert tests.test_rpacket_get_virtual_packet_flag(int_value)

@pytest.mark.parametrize(
	"int_value", integers
	)
def test_rpacket_get_virtual_packet(int_value):
	assert tests.test_rpacket_get_virtual_packet(int_value)

@pytest.mark.parametrize(
	"int_value", integers
	)
def test_rpacket_get_next_shell_id(int_value):
	assert tests.test_rpacket_get_next_shell_id(int_value)


# Testing functions without any parameters
def test_rpacket_get_last_line():
	assert tests.test_rpacket_get_last_line()

def test_rpacket_get_close_line():
	assert tests.test_rpacket_get_close_line()

def test_rpacket_get_status():
	assert tests.test_rpacket_get_status()