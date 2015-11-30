import os
import random
from ctypes import CDLL, c_double

import pytest
import numpy as np

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

pytestmark = pytest.mark.skipif(True, reason='problem with the files')

tests = CDLL(test_path)

np.random.seed(1)


def get_doubles(size=10):
    from sys import float_info as floats
    MAX_FLOAT, MIN_FLOAT = floats[0], floats[3]
    return map(c_double, np.random.uniform(
        MIN_FLOAT, MAX_FLOAT, size=10))

@pytest.fixture(params=get_doubles())
def double_value(request):
    return request.param


def get_integers(size=10):
    MAX_INT, MIN_INT = 10000, -10000
    return np.random.randint(
        MIN_INT, MAX_INT, size=10)

@pytest.fixture(params=get_integers())
def int_value(request):
    return request.param


def get_unsigned_integers(size=10):
    MAX_INT = 10000
    return np.random.randint(
        0, MAX_INT, size=10)

@pytest.fixture(params=get_unsigned_integers())
def unsigned_int_value(request):
    return request.param


# Testing functions with float(C double) valued parameters

def test_rpacket_get_nu(double_value):
    assert tests.test_rpacket_get_nu(double_value)


def test_rpacket_get_mu(double_value):
    assert tests.test_rpacket_get_mu(double_value)


def test_rpacket_get_energy(double_value):
    assert tests.test_rpacket_get_energy(double_value)


def test_rpacket_get_r(double_value):
    assert tests.test_rpacket_get_r(double_value)


def test_rpacket_get_tau_event(double_value):
    assert tests.test_rpacket_get_tau_event(double_value)


def test_rpacket_get_nu_line(double_value):
    assert tests.test_rpacket_get_nu_line(double_value)


def test_rpacket_get_d_boundary(double_value):
    assert tests.test_rpacket_get_d_boundary(double_value)


def test_rpacket_get_d_electron(double_value):
    assert tests.test_rpacket_get_d_electron(double_value)


def test_rpacket_get_d_line(double_value):
    assert tests.test_rpacket_get_d_line(double_value)


# Testing functions with Unsigned Integer valued parameters

def test_rpacket_get_current_shell_id(unsigned_int_value):
    assert tests.test_rpacket_get_current_shell_id(unsigned_int_value)


def test_rpacket_get_next_line_id(unsigned_int_value):
    assert tests.test_rpacket_get_next_line_id(unsigned_int_value)


# Testing functions with Integer valued parameters

def test_rpacket_get_recently_crossed_boundary(int_value):
    assert tests.test_rpacket_get_recently_crossed_boundary(int_value)


def test_rpacket_get_virtual_packet_flag(int_value):
    assert tests.test_rpacket_get_virtual_packet_flag(int_value)


def test_rpacket_get_virtual_packet(int_value):
    assert tests.test_rpacket_get_virtual_packet(int_value)


def test_rpacket_get_next_shell_id(int_value):
    assert tests.test_rpacket_get_next_shell_id(int_value)


# Testing functions without any parameters
def test_rpacket_get_last_line():
    assert tests.test_rpacket_get_last_line()

def test_rpacket_get_close_line():
    assert tests.test_rpacket_get_close_line()

def test_rpacket_get_status():
    assert tests.test_rpacket_get_status()
