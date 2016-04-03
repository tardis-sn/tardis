import os
import pytest
from ctypes import CDLL, byref, c_uint, c_int64, c_double, c_ulong
from numpy.testing import assert_almost_equal

from tardis import __path__ as path
from tardis.montecarlo.struct import RPacket, StorageModel, RKState

test_path = os.path.join(path[0], 'montecarlo', 'test_montecarlo.so')
cmontecarlo_tests = CDLL(test_path)

cmontecarlo_filepath = os.path.join(path[0], 'montecarlo', 'montecarlo.so')
cmontecarlo_methods = CDLL(cmontecarlo_filepath)


@pytest.fixture(scope="function")
def packet():
    packet_default = {
        'nu': 0.4,
        'mu': 0.3,
        'energy': 0.9,
        'r': 7.5e14,
        'tau_event': 2.9e13,
        'nu_line': 0.2,
        'current_shell_id': 0,
        'next_line_id': 1,
        'last_line': 0,
        'close_line': 0,
        'current_continuum_id': 1,
        'virtual_packet_flag': 1,
        'virtual_packet': 0,
        'status': 0,
        'id': 0
    }
    return RPacket(**packet_default)


@pytest.fixture(scope="function")
def model():
    model_default = {
        'last_line_interaction_in_id': (c_int64 * 2)(*([0] * 2)),
        'last_line_interaction_shell_id': (c_int64 * 2)(*([0] * 2)),
        'last_line_interaction_type': (c_int64 * 2)(*([2])),

        'no_of_shells': 2,

        'r_inner': (c_double * 2)(*[6.912e14, 8.64e14]),
        'r_outer': (c_double * 2)(*[8.64e14, 1.0368e15]),

        'time_explosion': 5.2e7,
        'inverse_time_explosion': 1 / 5.2e7,

        'electron_densities': (c_double * 2)(*[1.0e9] * 2),
        'inverse_electron_densities': (c_double * 2)(*[1.0e-9] * 2),

        'line_list_nu': (c_double * 5)(*[1.26318289e+16, 1.26318289e+16,
                                         1.23357675e+16, 1.23357675e+16,
                                         1.16961598e+16]),

        'continuum_list_nu': (c_double * 20000)(*([1.e13] * 20000)),

        'line_lists_tau_sobolevs': (c_double * 1000)(*([1.e-5] * 1000)),
        'line_lists_j_blues': (c_double * 2)(*([1.e-10] * 2)),
        'line_lists_j_blues_nd': 0,

        'no_of_lines': 2,
        'no_of_edges': 100,

        'line_interaction_id': 0,
        'line2macro_level_upper': (c_int64 * 2)(*([0] * 2)),

        'js': (c_double * 2)(*([0.0] * 2)),
        'nubars': (c_double * 2)(*([0.0] * 2)),

        'spectrum_start_nu': 1.e14,
        'spectrum_delta_nu': 293796608840.0,
        'spectrum_end_nu': 6.e15,

        'spectrum_virt_start_nu': 1e14,
        'spectrum_virt_end_nu': 6e15,
        'spectrum_virt_nu': (c_double * 20000)(*([0.0] * 20000)),

        'sigma_thomson': 6.652486e-25,
        'inverse_sigma_thomson': 1 / 6.652486e-25,

        'inner_boundary_albedo': 0.0,
        'reflective_inner_boundary': 0,

        'chi_bf_tmp_partial': (c_double * 20000)(*([160.0] * 20000)),
        't_electrons': (c_double * 2)(*([0.0] * 2)),

        'l_pop': (c_double * 20000)(*([2.0] * 20000)),
        'l_pop_r': (c_double * 20000)(*([3.0] * 20000)),
    }
    return StorageModel(**model_default)


@pytest.fixture(scope="function")
def mt_state():
    mt_state_default = {
        'key': (c_ulong * 624)(*([0] * 624)),
        'pos': 0,
        'has_gauss': 0,
        'gauss': 0.0
    }
    return RKState(**mt_state_default)


@pytest.mark.parametrize(
    ['packet_params', 'model_params', 'expected'],
    [({'mu': 0.3, 'r': 7.5e14}, {'inverse_time_explosion': 1 / 5.2e7},
      0.9998556693818854),

     ({'mu': -.3, 'r': 8.1e14}, {'inverse_time_explosion': 1 / 2.6e7},
      1.0003117541351274)]
)
def test_rpacket_doppler_factor(packet_params, model_params, expected, packet, model):
    packet.mu = packet_params['mu']
    packet.r = packet_params['r']
    model.inverse_time_explosion = model_params['inverse_time_explosion']

    cmontecarlo_methods.rpacket_doppler_factor.restype = c_double
    obtained = cmontecarlo_methods.rpacket_doppler_factor(byref(packet), byref(model))

    assert_almost_equal(obtained, expected)


@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    [({'mu': 0.3, 'r': 7.5e14},
      {'d_boundary': 259376919351035.88}),

     ({'mu': -.3, 'r': 7.5e13},
      {'d_boundary': -664987228972291.5}),

     ({'mu': -.3, 'r': 7.5e14},
      {'d_boundary': 709376919351035.9})]
)
def test_compute_distance2boundary(packet_params, expected_params, packet, model):
    packet.mu = packet_params['mu']
    packet.r = packet_params['r']

    cmontecarlo_methods.compute_distance2boundary(byref(packet), byref(model))

    assert_almost_equal(packet.d_boundary, expected_params['d_boundary'])


@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    [({'nu_line': 0.1, 'next_line_id': 0, 'last_line': 1},
      {'tardis_error': 0, 'd_line': 1e+99}),

     ({'nu_line': 0.2, 'next_line_id': 1, 'last_line': 0},
      {'tardis_error': 0, 'd_line': 7.792353908000001e+17}),

     ({'nu_line': 0.5, 'next_line_id': 1, 'last_line': 0},
      {'tardis_error': 2, 'd_line': 0.0}),

     ({'nu_line': 0.6, 'next_line_id': 0, 'last_line': 0},
      {'tardis_error': 2, 'd_line': 0.0})]
)
def test_compute_distance2line(packet_params, expected_params, packet, model):
    packet.nu_line = packet_params['nu_line']
    packet.next_line_id = packet_params['next_line_id']
    packet.last_line = packet_params['last_line']

    packet.d_line = 0.0
    cmontecarlo_methods.compute_distance2line.restype = c_uint
    obtained_tardis_error = cmontecarlo_methods.compute_distance2line(byref(packet), byref(model))

    assert_almost_equal(packet.d_line, expected_params['d_line'])
    assert obtained_tardis_error == expected_params['tardis_error']


@pytest.mark.parametrize(
    ['packet_params', 'model_params', 'expected_params'],
    [({'virtual_packet': 0},
      {'cont_status': 0},
      {'chi_cont': 6.652486e-16, 'd_cont': 4.359272608766106e+28}),

     ({'virtual_packet': 1},
      {'cont_status': 0},
      {'chi_cont': 6.652486e-16, 'd_cont': 1e+99})]
)
def test_compute_distance2continuum(packet_params, model_params,
                                    expected_params, packet, model):
    packet.virtual_packet = packet_params['virtual_packet']
    model.cont_status = model_params['cont_status']

    cmontecarlo_methods.compute_distance2continuum(byref(packet), byref(model))

    assert_almost_equal(packet.chi_cont, expected_params['chi_cont'])
    assert_almost_equal(packet.d_cont, expected_params['d_cont'])


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
