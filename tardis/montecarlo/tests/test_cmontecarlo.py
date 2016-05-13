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
import numpy as np
import pytest
from astropy import constants
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

# Wrap the shared object containing C methods, which are tested here.
cmontecarlo_filepath = os.path.join(path[0], 'montecarlo', 'montecarlo.so')
cmontecarlo_methods = CDLL(cmontecarlo_filepath)

# Constant for machine precision
EPS = np.finfo(np.float64).eps
# Constants to pertubate by machine precision
EPSP = 1 + EPS
EPSM = 1 - EPS

INVERSE_C = 1/constants.c.cgs.value


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
        next_shell_id=1,
        status=TARDIS_PACKET_STATUS_IN_PROCESS,
        id=0,
        chi_cont=6.652486e-16
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


"""
Important Tests:
----------------
The tests written further (till next block comment is encountered) have been
categorized as important tests, these tests correspond to methods which are
relatively old and stable code.
"""


@pytest.mark.parametrize(
    ['x', 'x_insert', 'imin', 'imax', 'expected_params'],
    [([5.0, 4.0, 3.0, 1.0], 2.0, 0, 3,
      {'result': 2, 'ret_val': TARDIS_ERROR_OK}),

     ([5.0, 4.0, 3.0, 2.0], 0.0, 0, 3,
      {'result': 0, 'ret_val': TARDIS_ERROR_BOUNDS_ERROR})]
)
def test_reverse_binary_search(x, x_insert, imin, imax, expected_params):
    x = (c_double * (imax - imin + 1))(*x)
    x_insert = c_double(x_insert)
    imin = c_int64(imin)
    imax = c_int64(imax)
    obtained_result = c_int64(0)

    cmontecarlo_methods.reverse_binary_search.restype = c_uint
    obtained_tardis_error = cmontecarlo_methods.reverse_binary_search(
                        byref(x), x_insert, imin, imax, byref(obtained_result))

    assert obtained_result.value == expected_params['result']
    assert obtained_tardis_error == expected_params['ret_val']


@pytest.mark.parametrize(
    ['nu', 'nu_insert', 'number_of_lines', 'expected_params'],
    [([0.5, 0.4, 0.3, 0.1], 0.2, 4,
      {'result': 3, 'ret_val': TARDIS_ERROR_OK}),

     ([0.5, 0.4, 0.3, 0.2], 0.1, 4,
      {'result': 4, 'ret_val': TARDIS_ERROR_OK}),

     ([0.4, 0.3, 0.2, 0.1], 0.5, 4,
      {'result': 0, 'ret_val': TARDIS_ERROR_OK})]
)
def test_line_search(nu, nu_insert, number_of_lines, expected_params):
    nu = (c_double * number_of_lines)(*nu)
    nu_insert = c_double(nu_insert)
    number_of_lines = c_int64(number_of_lines)
    obtained_result = c_int64(0)

    cmontecarlo_methods.line_search.restype = c_uint
    obtained_tardis_error = cmontecarlo_methods.line_search(
                        byref(nu), nu_insert, number_of_lines, byref(obtained_result))

    assert obtained_result.value == expected_params['result']
    assert obtained_tardis_error == expected_params['ret_val']


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
      {'tardis_error': TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE, 'd_line': np.nan}),

     ({'nu_line': 0.6, 'next_line_id': 0, 'last_line': 0},
      {'tardis_error': TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE, 'd_line': np.nan})]
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
    ['packet_params', 'expected_params'],
    [({'nu': 0.4, 'mu': 0.3, 'energy': 0.9, 'r': 7.5e14},
      {'mu': 0.3120599529139568, 'r': 753060422542573.9,
       'j': 8998701024436.969, 'nubar': 3598960894542.354}),

     ({'nu': 0.6, 'mu': -.5, 'energy': 0.5, 'r': 8.1e14},
      {'mu': -.4906548373534084, 'r': 805046582503149.2,
       'j': 5001298975563.031, 'nubar': 3001558973156.1387})]
)
def test_move_packet(packet_params, expected_params, packet, model):
    packet.nu = packet_params['nu']
    packet.mu = packet_params['mu']
    packet.energy = packet_params['energy']
    packet.r = packet_params['r']

    cmontecarlo_methods.move_packet(byref(packet), byref(model), c_double(1.e13))

    assert_almost_equal(packet.mu, expected_params['mu'])
    assert_almost_equal(packet.r, expected_params['r'])

    assert_almost_equal(model.js[packet.current_shell_id], expected_params['j'])
    assert_almost_equal(model.nubars[packet.current_shell_id], expected_params['nubar'])


# This function can be parametrized using hypothesis now
# It should be easy to add new parameters with random values
@pytest.mark.parametrize(
    ['distance', 'j_blue_idx'],
    [(1e13, 0),
     (2e13, 0),
     (5e13, 1),
     (1e14, 1)]
)
def test_increment_j_blue_estimator(distance, j_blue_idx, packet, model):
    r_interaction = np.sqrt(
            packet.r**2 + distance**2 +
            2.0 * packet.r * distance * packet.mu)
    mu_interaction = (packet.mu * packet.r + distance) / r_interaction
    doppler_factor = (
            1.0 - mu_interaction * r_interaction *
            model.inverse_time_explosion * INVERSE_C)
    comov_energy = packet.energy * doppler_factor

    model.line_lists_j_blues[j_blue_idx] = 0
    cmontecarlo_methods.increment_j_blue_estimator(
            byref(packet), byref(model),
            c_double(distance), c_int64(j_blue_idx))

    assert_almost_equal(
            model.line_lists_j_blues[j_blue_idx],
            comov_energy / packet.nu)


@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    [({'virtual_packet': 0, 'current_shell_id': 0, 'next_shell_id': 1},
      {'status': TARDIS_PACKET_STATUS_IN_PROCESS, 'current_shell_id': 1}),

     ({'virtual_packet': 1, 'current_shell_id': 1, 'next_shell_id': 1},
      {'status': TARDIS_PACKET_STATUS_EMITTED, 'current_shell_id': 1,
       'tau_event': 29000000000000.008}),

     ({'virtual_packet': 1, 'current_shell_id': 0, 'next_shell_id': -1},
      {'status': TARDIS_PACKET_STATUS_REABSORBED, 'current_shell_id': 0,
       'tau_event': 29000000000000.008})]
)
def test_move_packet_across_shell_boundary(packet_params, expected_params,
                                           packet, model, mt_state):
    packet.virtual_packet = packet_params['virtual_packet']
    packet.current_shell_id = packet_params['current_shell_id']
    packet.next_shell_id = packet_params['next_shell_id']

    cmontecarlo_methods.move_packet_across_shell_boundary(byref(packet), byref(model),
                                                          c_double(1.e13), byref(mt_state))

    if packet_params['virtual_packet'] == 1:
        assert_almost_equal(packet.tau_event, expected_params['tau_event'])
    assert packet.status == expected_params['status']
    assert packet.current_shell_id == expected_params['current_shell_id']


@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    [({'nu': 0.4, 'mu': 0.3, 'energy': 0.9, 'r': 7.5e14},
      {'nu': 0.39974659819356556, 'energy': 0.8994298459355226}),

     ({'nu': 0.6, 'mu': -.5, 'energy': 0.5, 'r': 8.1e14},
      {'nu': 0.5998422620533325, 'energy': 0.4998685517111104})]
)
def test_montecarlo_thomson_scatter(packet_params, expected_params, packet,
                                   model, mt_state):
    packet.nu = packet_params['nu']
    packet.mu = packet_params['mu']
    packet.energy = packet_params['energy']
    packet.r = packet_params['r']

    cmontecarlo_methods.montecarlo_thomson_scatter(byref(packet), byref(model),
                                                   c_double(1.e13), byref(mt_state))

    assert_almost_equal(packet.nu, expected_params['nu'])
    assert_almost_equal(packet.energy, expected_params['energy'])


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



"""
Difficult Tests:
----------------
The tests written further are more complex than previous tests. They require
proper design procedure. They are not taken up yet and intended to be
completed together in future.
"""


@pytest.mark.skipif(True, reason="Yet to be written.")
def test_montecarlo_one_packet(packet, model, mt_state):
    pass


@pytest.mark.skipif(True, reason="Yet to be written.")
def test_montecarlo_one_packet_loop(packet, model, mt_state):
    pass


@pytest.mark.skipif(True, reason="Yet to be written.")
def test_montecarlo_main_loop(packet, model, mt_state):
    pass


@pytest.mark.skipif(True, reason="Yet to be written.")
def test_montecarlo_event_handler(packet, model, mt_state):
    pass


"""
Not Yet Relevant Tests:
-----------------------
The tests written further (till next block comment is encountered) are for the
methods related to Continuum interactions. These are not required to be tested
on current master and can be skipped for now.
"""


@pytest.mark.skipif(True, reason="Not yet relevant")
@pytest.mark.parametrize(
    ['packet_params', 'expected'],
    [({'nu': 0.1, 'mu': 0.3, 'r': 7.5e14}, 2.5010827921809502e+26),
     ({'nu': 0.2, 'mu': -.3, 'r': 7.7e14}, 3.123611229395459e+25)]
)
def test_bf_cross_section(packet_params, expected, packet, model):
    packet.nu = packet_params['nu']
    packet.mu = packet_params['mu']
    packet.r = packet_params['r']

    cmontecarlo_methods.rpacket_doppler_factor.restype = c_double
    doppler_factor = cmontecarlo_methods.rpacket_doppler_factor(byref(packet), byref(model))
    comov_nu = packet.nu * doppler_factor

    cmontecarlo_methods.bf_cross_section.restype = c_double
    obtained = cmontecarlo_methods.bf_cross_section(byref(model), c_int64(0),
                                                    c_double(comov_nu))

    assert_almost_equal(obtained, expected)


# TODO: fix underlying method and update expected values in testcases.
# For loop is not being executed in original method, and hence bf_helper
# always remains zero. Reason for for loop not executed:
#         "current_continuum_id = no_of_continuum edges"
@pytest.mark.skipif(True, reason="Not yet relevant")
@pytest.mark.parametrize(
    ['packet_params', 'expected'],
    [({'nu': 0.1, 'mu': 0.3, 'r': 7.5e14}, 0.0),
     ({'nu': 0.2, 'mu': -.3, 'r': 7.7e14}, 0.0)]
)
def test_calculate_chi_bf(packet_params, expected, packet, model):
    packet.nu = packet_params['nu']
    packet.mu = packet_params['mu']
    packet.r = packet_params['r']

    cmontecarlo_methods.calculate_chi_bf(byref(packet), byref(model))

    assert_almost_equal(packet.chi_bf, expected)


@pytest.mark.skipif(True, reason="Not yet relevant")
def test_montecarlo_continuum_event_handler(packet_params, continuum_status, expected,
                                            packet, model, mt_state):
    packet.chi_cont = packet_params['chi_cont']
    packet.chi_th = packet_params['chi_th']
    packet.chi_bf = packet.chi_cont - packet.chi_th
    model.cont_status = continuum_status

    obtained = cmontecarlo_methods.montecarlo_continuum_event_handler(byref(packet),
                                                      byref(model), byref(mt_state))


@pytest.mark.skipif(True, reason="Not yet relevant")
def test_montecarlo_free_free_scatter(packet, model, mt_state):
    cmontecarlo_methods.montecarlo_free_free_scatter(byref(packet), byref(model),
                                                     c_double(1.e13), byref(mt_state))

    assert_equal(packet.status, TARDIS_PACKET_STATUS_REABSORBED)


@pytest.mark.skipif(True, reason="Not yet relevant")
def test_montecarlo_bound_free_scatter(packet, model, mt_state):
    cmontecarlo_methods.montecarlo_bound_free_scatter(byref(packet), byref(model),
                                                     c_double(1.e13), byref(mt_state))

    assert_equal(packet.status, TARDIS_PACKET_STATUS_REABSORBED)
