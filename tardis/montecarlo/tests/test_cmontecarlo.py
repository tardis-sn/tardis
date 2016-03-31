import os
from ctypes import CDLL, c_double
import pytest
from numpy.testing import assert_almost_equal

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'test_montecarlo.so')
cmontecarlo_tests = CDLL(test_path)


def test_compute_distance2boundary():
    distance_to_boundary = 259376919351035.88
    cmontecarlo_tests.test_compute_distance2boundary.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_compute_distance2boundary(),
                        distance_to_boundary)


def test_compute_distance2line():
    distance_to_line = 7.792353908000001e+17
    cmontecarlo_tests.test_compute_distance2line.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_compute_distance2line(),
                        distance_to_line)


def test_compute_distance2continuum():
    distance_to_electron = 4.359272608766106e+28
    cmontecarlo_tests.test_compute_distance2continuum.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_compute_distance2continuum(),
                        distance_to_electron)


def test_rpacket_doppler_factor():
    doppler_factor = 0.9998556693818854
    cmontecarlo_tests.test_rpacket_doppler_factor.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_rpacket_doppler_factor(),
                        doppler_factor)


@pytest.mark.skipif(True, reason='Bad test design')
def test_move_packet():
    doppler_factor = 0.9998556693818854
    cmontecarlo_tests.test_move_packet.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_move_packet(),
                        doppler_factor)


def test_increment_j_blue_estimator():
    j_blue = 1.1249855669381885
    cmontecarlo_tests.test_increment_j_blue_estimator.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_increment_j_blue_estimator(),
                        j_blue)


def test_montecarlo_line_scatter():
    assert cmontecarlo_tests.test_montecarlo_line_scatter()


def test_move_packet_across_shell_boundary():
    assert cmontecarlo_tests.test_move_packet_across_shell_boundary()


def test_montecarlo_one_packet():
    assert cmontecarlo_tests.test_montecarlo_one_packet()


def test_montecarlo_one_packet_loop():
    assert cmontecarlo_tests.test_montecarlo_one_packet_loop() == 0


def test_montecarlo_thomson_scatter():
    assert cmontecarlo_tests.test_montecarlo_thomson_scatter()


def test_calculate_chi_bf():
    chi_bf = 1.0006697327643788
    cmontecarlo_tests.test_calculate_chi_bf.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_calculate_chi_bf(),
                        chi_bf)


@pytest.mark.xfail
def test_montecarlo_bound_free_scatter():
    assert cmontecarlo_tests.test_montecarlo_bound_free_scatter() == 1


@pytest.mark.xfail
def test_bf_cross_section():
    bf_cross_section = 0.0
    cmontecarlo_tests.test_bf_cross_section.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_bf_cross_section(),
                        bf_cross_section)


def test_montecarlo_free_free_scatter():
    assert cmontecarlo_tests.test_montecarlo_free_free_scatter() == 2
