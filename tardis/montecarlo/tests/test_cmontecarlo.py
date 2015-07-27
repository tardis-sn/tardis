import os
import random
from ctypes import CDLL, c_double

import pytest
import numpy as np
from numpy.testing import assert_almost_equal

from tardis import __path__ as path

test_path = os.path.join(path[0], 'montecarlo', 'montecarlo.so')

tests = CDLL(test_path)

tests.init_rpacket()
tests.init_storage_model()

def test_compute_distance2boundary():
	distance_to_boundary = 259376919351035.88
	tests.test_compute_distance2boundary.restype = c_double
	assert_almost_equal(tests.test_compute_distance2boundary(), 
		distance_to_boundary)

def test_compute_distance2line():
	distance_to_line = 7.792353908000001e+17
	tests.test_compute_distance2line.restype = c_double
	assert_almost_equal(tests.test_compute_distance2line(), 
		distance_to_line)

def test_compute_distance2continuum():
	distance_to_electron = 4.359272608766106e+28
	tests.test_compute_distance2continuum.restype = c_double
	assert_almost_equal(tests.test_compute_distance2continuum(), 
		distance_to_electron)

def test_rpacket_doppler_factor():
	doppler_factor = 0.9998556693818854
	tests.test_rpacket_doppler_factor.restype = c_double
	assert_almost_equal(tests.test_rpacket_doppler_factor(), 
		doppler_factor)

def test_move_packet():
	doppler_factor = 0.9998556693818854
	tests.test_move_packet.restype = c_double
	assert_almost_equal(tests.test_move_packet(),
		doppler_factor)

def test_increment_j_blue_estimator():
	j_blue = 1.1249855669381885
	tests.test_increment_j_blue_estimator.restype = c_double
	assert_almost_equal( tests.test_increment_j_blue_estimator(),
		j_blue)

def test_montecarlo_line_scatter():
	assert tests.test_montecarlo_line_scatter()

def test_move_packet_across_shell_boundary():
	assert tests.test_move_packet_across_shell_boundary()

def test_montecarlo_one_packet():
	assert tests.test_montecarlo_one_packet()

def test_montecarlo_one_packet_loop():
	assert tests.test_montecarlo_one_packet_loop() == 0

def test_montecarlo_thomson_scatter():
	assert tests.test_montecarlo_thomson_scatter()

def test_calculate_chi_bf():
	chi_bf = 1.0006697327643788
	tests.test_calculate_chi_bf.restype = c_double
	assert_almost_equal(tests.test_calculate_chi_bf(),
		chi_bf)

def test_montecarlo_bound_free_scatter():
	assert tests.test_montecarlo_bound_free_scatter() == 1

def test_bf_cross_section():
	bf_cross_section = -1
	tests.test_bf_cross_section.restype = c_double
	assert_almost_equal(tests.test_bf_cross_section(),
		bf_cross_section)

def test_montecarlo_free_free_scatter():
	assert tests.test_montecarlo_free_free_scatter() == 2

# def test_macro_atom():
#	assert tests.test_macro_atom()