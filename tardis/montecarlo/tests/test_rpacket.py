import os
import random
import inspect
from ctypes import CDLL, c_double

import pytest
import numpy as np

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

tests = CDLL(test_path)

SEED = 26061963

def get_random_float(size=10):
    from sys import float_info as floats
    MAX_FLOAT, MIN_FLOAT = floats[0], floats[3]
    np.random.seed(SEED)
    return np.random.uniform(MIN_FLOAT, MAX_FLOAT, size=size)

@pytest.fixture(params=get_random_float())
def random_double_value(request):
    return c_double(request.param)


def get_random_unsigned_int(size=10):
    MAX_INT, MIN_INT = 10000, -10000
    np.random.seed(SEED)
    return np.random.randint(0, MAX_INT, size=size)

@pytest.fixture(params=get_random_unsigned_int())
def random_unsigned_int(request):
    return request.param

def get_random_int(size=10):
    MAX_INT, MIN_INT = 10000, -10000
    np.random.seed(SEED)
    return np.random.randint(MIN_INT, MAX_INT, size=size)

@pytest.fixture(params=get_random_int())
def random_int(request):
    return request.param




# Testing functions with float(C double) valued parameters

def test_rpacket_get_nu(random_double_value):
    assert tests.test_rpacket_get_nu(random_double_value)


def test_rpacket_get_mu(random_double_value):
    assert tests.test_rpacket_get_mu(random_double_value)


def test_rpacket_get_energy(random_double_value):
    assert tests.test_rpacket_get_energy(random_double_value)


def test_rpacket_get_r(random_double_value):
    assert tests.test_rpacket_get_r(random_double_value)


def test_rpacket_get_tau_event(random_double_value):
    assert tests.test_rpacket_get_tau_event(random_double_value)


def test_rpacket_get_nu_line(random_double_value):
    assert tests.test_rpacket_get_nu_line(random_double_value)


def test_rpacket_get_d_boundary(random_double_value):
    assert tests.test_rpacket_get_d_boundary(random_double_value)


def test_rpacket_get_d_electron(random_double_value):
    assert tests.test_rpacket_get_d_electron(random_double_value)

def test_rpacket_get_d_line(random_double_value):
    assert tests.test_rpacket_get_d_line(random_double_value)


# Testing functions with Unsigned Integer valued parameters

def test_rpacket_get_current_shell_id(random_unsigned_int):
    assert tests.test_rpacket_get_current_shell_id(random_unsigned_int)


def test_rpacket_get_next_line_id(random_unsigned_int):
    assert tests.test_rpacket_get_next_line_id(random_unsigned_int)


# Testing functions with Integer valued parameters
def test_rpacket_get_recently_crossed_boundary(random_int):
    assert tests.test_rpacket_get_recently_crossed_boundary(random_int)

def test_rpacket_get_virtual_packet_flag(random_int):
    assert tests.test_rpacket_get_virtual_packet_flag(random_int)

def test_rpacket_get_virtual_packet(random_int):
    assert tests.test_rpacket_get_virtual_packet(random_int)

def test_rpacket_get_next_shell_id(random_int):
    assert tests.test_rpacket_get_next_shell_id(random_int)


# Testing functions without any parameters
def test_rpacket_get_last_line():
    assert tests.test_rpacket_get_last_line()

def test_rpacket_get_close_line():
    assert tests.test_rpacket_get_close_line()

def test_rpacket_get_status():
    assert tests.test_rpacket_get_status()