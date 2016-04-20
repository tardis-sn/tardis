"""
Unit tests for methods in `tardis/montecarlo/src/cmontecarlo.c`.
* `ctypes` library is used to wrap C methods and expose them to python.


Probable Reasons for Failing Tests:
-----------------------------------

1. Change made in C struct declarations:
  - Reflect the changes done in C structs, into Python counterparts.
  - Check **tardis/montecarlo/struct.py**.

2. Return type of any method changed:
  - Modify the `restype` parameter in the test method here.
  - For example:
        ```
        cmontecarlo_methods.rpacket_doppler_factor.restype = c_double
        ```

3. Underlying logic modified:
  - Check whether the changes made in C method are logically correct.
  - If the changes made were correct and necessary, update the corresponding
    test case.


General Test Design Procedure:
------------------------------

Please follow this design procedure while adding a new test:

1. Parametrization as per desire of code coverage.
  - C tests have different flows controlled by conditional statements.
    Parameters checked in conditions can be provided in different testcases.
  - Keep consistency with variable names as (in order):
    - `packet_params`
    - `model_params`
    - `expected_params` (`expected` if only one value to be asserted.)
  - Suggested variable names can be compromised if readability of the test
    increases.

2. Test Method body:
  - Keep name as `test_` + `(name of C method)`.
  - Refer to method `test_rpacket_doppler_factor` below for description.
"""

import os
import pytest
from ctypes import CDLL, byref, c_uint, c_int64, c_double, c_ulong
from numpy.testing import assert_equal, assert_almost_equal

from tardis import __path__ as path
from tardis.montecarlo.struct import (
    RPacket, StorageModel, RKState,
    TARDIS_ERROR_OK,
    TARDIS_ERROR_BOUNDS_ERROR,
    TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE,
    TARDIS_PACKET_STATUS_IN_PROCESS,
    TARDIS_PACKET_STATUS_EMITTED,
    TARDIS_PACKET_STATUS_REABSORBED,
    CONTINUUM_OFF,
    CONTINUUM_ON
)

# Wrap the shared object containing tests for C methods, written in C.
# TODO: Shift all tests here in Python and completely remove this test design.
test_path = os.path.join(path[0], 'montecarlo', 'test_montecarlo.so')
cmontecarlo_tests = CDLL(test_path)

# Wrap the shared object containing C methods, which are tested here.
cmontecarlo_filepath = os.path.join(path[0], 'montecarlo', 'montecarlo.so')
cmontecarlo_methods = CDLL(cmontecarlo_filepath)


@pytest.fixture(scope="function")
def packet():
    """Fixture to return `RPacket` object with default params initialized."""
    return RPacket(
        nu=0.4,
        mu=0.3,
        energy=0.9,
        r=7.5e14,
        tau_event=2.9e13,
        nu_line=0.2,
        current_shell_id=0,
        next_line_id=1,
        last_line=0,
        close_line=0,
        current_continuum_id=1,
        virtual_packet_flag=1,
        virtual_packet=0,
        status=TARDIS_PACKET_STATUS_IN_PROCESS,
        id=0
    )


@pytest.fixture(scope="function")
def model():
    """Fixture to return `StorageModel` object with default params initialized."""
    return StorageModel(
        last_line_interaction_in_id=(c_int64 * 2)(*([0] * 2)),
        last_line_interaction_shell_id=(c_int64 * 2)(*([0] * 2)),
        last_line_interaction_type=(c_int64 * 2)(*([2])),

        no_of_shells=2,

        r_inner=(c_double * 2)(*[6.912e14, 8.64e14]),
        r_outer=(c_double * 2)(*[8.64e14, 1.0368e15]),

        time_explosion=5.2e7,
        inverse_time_explosion=1 / 5.2e7,

        electron_densities=(c_double * 2)(*[1.0e9] * 2),
        inverse_electron_densities=(c_double * 2)(*[1.0e-9] * 2),

        line_list_nu=(c_double * 5)(*[1.26318289e+16, 1.26318289e+16,
                                         1.23357675e+16, 1.23357675e+16,
                                         1.16961598e+16]),

        continuum_list_nu=(c_double * 20000)(*([1.e13] * 20000)),

        line_lists_tau_sobolevs=(c_double * 1000)(*([1.e-5] * 1000)),
        line_lists_j_blues=(c_double * 2)(*([1.e-10] * 2)),
        line_lists_j_blues_nd=0,

        no_of_lines=2,
        no_of_edges=100,

        line_interaction_id=0,
        line2macro_level_upper=(c_int64 * 2)(*([0] * 2)),

        js=(c_double * 2)(*([0.0] * 2)),
        nubars=(c_double * 2)(*([0.0] * 2)),

        spectrum_start_nu=1.e14,
        spectrum_delta_nu=293796608840.0,
        spectrum_end_nu=6.e15,

        spectrum_virt_start_nu=1e14,
        spectrum_virt_end_nu=6e15,
        spectrum_virt_nu=(c_double * 20000)(*([0.0] * 20000)),

        sigma_thomson=6.652486e-25,
        inverse_sigma_thomson=1 / 6.652486e-25,

        inner_boundary_albedo=0.0,
        reflective_inner_boundary=0,

        chi_bf_tmp_partial=(c_double * 20000)(*([160.0] * 20000)),
        t_electrons=(c_double * 2)(*([0.0] * 2)),

        l_pop=(c_double * 20000)(*([2.0] * 20000)),
        l_pop_r=(c_double * 20000)(*([3.0] * 20000)),
        cont_status=CONTINUUM_OFF
    )


@pytest.fixture(scope="function")
def mt_state():
    """Fixture to return `RKState` object with default params initialized."""
    return RKState(
        key=(c_ulong * 624)(*([0] * 624)),
        pos=0,
        has_gauss=0,
        gauss=0.0
    )


@pytest.mark.parametrize(
    ['mu', 'r', 'inv_t_exp', 'expected'],
    [(0.3, 7.5e14, 1 / 5.2e7, 0.9998556693818854),
     (-.3, 8.1e14, 1 / 2.6e7, 1.0003117541351274)]
)
def test_rpacket_doppler_factor(mu, r, inv_t_exp, expected, packet, model):
    # Set the params from test cases here
    packet.mu = mu
    packet.r = r
    model.inverse_time_explosion = inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    # Set `restype` attribute if returned quantity is used
    cmontecarlo_methods.rpacket_doppler_factor.restype = c_double
    # Call the C method (make sure to pass quantities as `ctypes` data types)
    obtained = cmontecarlo_methods.rpacket_doppler_factor(byref(packet), byref(model))

    # Perform required assertions
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


# TODO: split this into two tests - one to assert errors and other for d_line
@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    [({'nu_line': 0.1, 'next_line_id': 0, 'last_line': 1},
      {'tardis_error': TARDIS_ERROR_OK, 'd_line': 1e+99}),

     ({'nu_line': 0.2, 'next_line_id': 1, 'last_line': 0},
      {'tardis_error': TARDIS_ERROR_OK, 'd_line': 7.792353908000001e+17}),

     ({'nu_line': 0.5, 'next_line_id': 1, 'last_line': 0},
      {'tardis_error': TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE, 'd_line': 0.0}),

     ({'nu_line': 0.6, 'next_line_id': 0, 'last_line': 0},
      {'tardis_error': TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE, 'd_line': 0.0})]
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
    ['packet_params', 'expected_params'],
    [({'virtual_packet': 0},
     {'chi_cont': 6.652486e-16, 'd_cont': 4.359272608766106e+28}),

     ({'virtual_packet': 1},
      {'chi_cont': 6.652486e-16, 'd_cont': 1e+99})]
)
def test_compute_distance2continuum(packet_params, expected_params, packet, model):
    packet.virtual_packet = packet_params['virtual_packet']

    cmontecarlo_methods.compute_distance2continuum(byref(packet), byref(model))

    assert_almost_equal(packet.chi_cont, expected_params['chi_cont'])
    assert_almost_equal(packet.d_cont, expected_params['d_cont'])


@pytest.mark.parametrize(
    ['packet_params', 'j_blue_idx', 'expected'],
    [({'nu': 0.1, 'mu': 0.3, 'r': 7.5e14}, 0, 8.998643292289723),
     ({'nu': 0.2, 'mu': -.3, 'r': 7.7e14}, 0, 4.499971133976377),
     ({'nu': 0.5, 'mu': 0.5, 'r': 7.9e14}, 1, 0.719988453650551),
     ({'nu': 0.6, 'mu': -.5, 'r': 8.1e14}, 1, 0.499990378058792)]
)
def test_increment_j_blue_estimator(packet_params, j_blue_idx, expected, packet, model):
    packet.nu = packet_params['nu']
    packet.mu = packet_params['mu']
    packet.r = packet_params['r']

    cmontecarlo_methods.compute_distance2line(byref(packet), byref(model))
    cmontecarlo_methods.move_packet(byref(packet), byref(model), c_double(1.e13))
    cmontecarlo_methods.increment_j_blue_estimator(byref(packet), byref(model),
                                 c_double(packet.d_line), c_int64(j_blue_idx))

    assert_almost_equal(model.line_lists_j_blues[j_blue_idx], expected)


@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    # TODO: Add scientifically sound test cases.
    [({'virtual_packet': 1, 'tau_event': 2.9e13, 'last_line': 0},
      {'tau_event': 2.9e13, 'next_line_id': 2}),

     ({'virtual_packet': 0, 'tau_event': 2.9e13, 'last_line': 0},
      {'tau_event': 2.9e13, 'next_line_id': 2}),

     ({'virtual_packet': 0, 'tau_event': 2.9e13, 'last_line': 0},
      {'tau_event': 2.9e13, 'next_line_id': 2}),
     ]
)
def test_montecarlo_line_scatter(packet_params, expected_params, packet, model, mt_state):
    packet.virtual_packet = packet_params['virtual_packet']
    packet.tau_event = packet_params['tau_event']
    packet.last_line = packet_params['last_line']

    cmontecarlo_methods.montecarlo_line_scatter(byref(packet), byref(model),
                                          c_double(1.e13), byref(mt_state))

    assert_almost_equal(packet.tau_event, expected_params['tau_event'])
    assert_almost_equal(packet.next_line_id, expected_params['next_line_id'])


def test_montecarlo_free_free_scatter(packet, model, mt_state):
    cmontecarlo_methods.montecarlo_free_free_scatter(byref(packet), byref(model),
                                                     c_double(1.e13), byref(mt_state))

    assert_equal(packet.status, TARDIS_PACKET_STATUS_REABSORBED)


# TODO: redesign subsequent tests according to tests written above.
@pytest.mark.skipif(True, reason='Bad test design')
def test_move_packet():
    doppler_factor = 0.9998556693818854
    cmontecarlo_tests.test_move_packet.restype = c_double
    assert_almost_equal(cmontecarlo_tests.test_move_packet(),
                        doppler_factor)


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
