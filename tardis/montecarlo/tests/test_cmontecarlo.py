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
import numpy as np
import pandas as pd
from ctypes import CDLL, byref, c_uint, c_int, c_int64, c_double, c_ulong, c_void_p, cast, POINTER, pointer, CFUNCTYPE
from numpy.testing import assert_equal, assert_almost_equal, assert_array_equal, assert_allclose

from tardis import __path__ as path
from tardis.montecarlo.struct import (
    RPacket, StorageModel, RKState,
    PhotoXsect1level,
    TARDIS_ERROR_OK,
    TARDIS_ERROR_BOUNDS_ERROR,
    TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE,
    TARDIS_PACKET_STATUS_IN_PROCESS,
    TARDIS_PACKET_STATUS_EMITTED,
    TARDIS_PACKET_STATUS_REABSORBED,
    CONTINUUM_OFF,
    CONTINUUM_ON,
    BF_EMISSION,
    FF_EMISSION,
    EXCITATION,
    IONIZATION,
    C, INVERSE_C
)

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
        next_shell_id=1,
        status=TARDIS_PACKET_STATUS_IN_PROCESS,
        id=0,
        chi_cont=6.652486e-16,
        chi_bf_tmp_partial=(c_double * 2)(),
        compute_chi_bf=True
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

        continuum_list_nu=(c_double * 2)(*([1.e13] * 2)),

        line_lists_tau_sobolevs=(c_double * 1000)(*([1.e-5] * 1000)),
        line_lists_j_blues=(c_double * 2)(*([1.e-10] * 2)),
        line_lists_j_blues_nd=0,

        line_lists_Edotlu=(c_double * 3)(*[0.0,0.0,1.0]), # Init to an explicit array

        no_of_lines=5,
        no_of_edges=2,

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

        chi_ff_factor=(c_double * 2)(*([1.0] * 2)),
        t_electrons=(c_double * 2)(*([1.0e4] * 2)),

        l_pop=(c_double * 20000)(*([2.0] * 20000)),
        l_pop_r=(c_double * 20000)(*([3.0] * 20000)),
        cont_status=CONTINUUM_OFF,
        ff_heating_estimator=(c_double * 2)(*([0.0] * 2)),
        fb_cooling_prob = (c_double * 2)(*([0.0] * 2)),
        ff_cooling_prob = (c_double *2)(*([0.0] * 2)),
        coll_exc_cooling_prob = (c_double *2)(*([0.0] * 2)),
        coll_ion_cooling_prob = (c_double *2)(*([0.0] * 2))
    )


@pytest.fixture(scope='module')
def continuum_compare_data_fname():
    fname = 'continuum_compare_data.hdf'
    return os.path.join(path[0], 'montecarlo', 'tests', 'data', fname)


@pytest.fixture(scope='module')
def continuum_compare_data(continuum_compare_data_fname, request):
    compare_data = pd.HDFStore(continuum_compare_data_fname, mode='r')

    def fin():
        compare_data.close()
    request.addfinalizer(fin)

    return compare_data


@pytest.fixture(scope="function")
def mt_state():
    """Fixture to return `RKState` object with default params initialized."""
    return RKState(
        key=(c_ulong * 624)(*([0] * 624)),
        pos=0,
        has_gauss=0,
        gauss=0.0
    )


@pytest.fixture(scope="function")
def mt_state_seeded(mt_state):
    seed = 23111963
    cmontecarlo_methods.rk_seed(seed, byref(mt_state))
    return mt_state


@pytest.fixture(scope="function")
def expected_ff_emissivity(continuum_compare_data):
    emissivities = continuum_compare_data['ff_emissivity']

    def ff_emissivity(t_electron):
        emissivity = emissivities[t_electron]
        nu_bins = emissivity['nu_bins'].values
        emissivity_value = emissivity['emissivity'].dropna().values

        return nu_bins, emissivity_value

    return ff_emissivity


@pytest.fixture(scope='module')
def get_rkstate(continuum_compare_data):
    data = continuum_compare_data['z2rkstate']

    def z2rkstate(z_random):
        key = (c_ulong * 624)(*data.loc[z_random, 'key'])
        pos = data.loc[z_random, 'pos']
        return RKState(
            key=key,
            pos=pos,
            has_gauss=0,
            gauss=0.0
        )

    return z2rkstate


@pytest.fixture(scope='function')
def model_w_edges(ion_edges, model):
    photo_xsect = (POINTER(PhotoXsect1level) * len(ion_edges))()

    for i, edge in enumerate(ion_edges):
        x_sect_1level = PhotoXsect1level()
        for key, value in edge.iteritems():
            if key in ['nu', 'x_sect']:
                value = (c_double * len(value))(*value)
            setattr(x_sect_1level, key, value)
        photo_xsect[i] = pointer(x_sect_1level)

    no_of_edges = len(ion_edges)
    continuum_list_nu = (c_double * no_of_edges)(*[edge['nu'][0] for edge in ion_edges])

    model.photo_xsect = photo_xsect
    model.continuum_list_nu = continuum_list_nu
    model.no_of_edges = no_of_edges

    estimator_size = model.no_of_shells * no_of_edges
    estims = ['photo_ion_estimator', 'stim_recomb_estimator', 'bf_heating_estimator', 'stim_recomb_cooling_estimator']
    for estimator in estims:
        setattr(model, estimator, (c_double * estimator_size)(*[0] * estimator_size))

    model.photo_ion_estimator_statistics = (c_int64 * estimator_size)(*[0] * estimator_size)
    return model


@pytest.fixture(scope='module')
def ion_edges():
    return [
        {'nu': [4.0e14, 4.1e14, 4.2e14, 4.3e14], 'x_sect': [1.0, 0.9, 0.8, 0.7], 'no_of_points': 4},
        {'nu': [3.0e14, 3.1e14, 3.2e14, 3.3e14, 3.4e14], 'x_sect': [1.0, 0.9, 0.8, 0.7, 0.6], 'no_of_points': 5},
        {'nu': [2.8e14, 3.0e14, 3.2e14, 3.4e14], 'x_sect': [2.0, 1.8, 1.6, 1.4], 'no_of_points': 4}
    ]


@pytest.fixture(scope='module')
def mock_sample_nu():
    SAMPLE_NUFUNC = CFUNCTYPE(c_double, POINTER(RPacket), POINTER(StorageModel), POINTER(RKState))

    def sample_nu_simple(packet, model, mt_state):
        return packet.contents.nu

    return SAMPLE_NUFUNC(sample_nu_simple)


def d_cont_setter(d_cont, model, packet):
    model.inverse_electron_densities[packet.current_shell_id] = c_double(1.0)
    model.inverse_sigma_thomson = c_double(1.0)
    packet.tau_event = c_double(d_cont)


def d_line_setter(d_line, model, packet):
    packet.mu = c_double(0.0)
    scale = d_line * 1e1
    model.time_explosion = c_double(INVERSE_C * scale)
    packet.nu = c_double(1.0)
    nu_line = (1. - d_line/scale)
    packet.nu_line = c_double(nu_line)


def d_boundary_setter(d_boundary, model, packet):
    packet.mu = c_double(1e-16)
    r_outer = 2. * d_boundary
    model.r_outer[packet.current_shell_id] = r_outer

    r = np.sqrt(r_outer**2 - d_boundary**2)
    packet.r = r


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


@pytest.mark.parametrize(
    ['distances', 'expected'],
    [({'boundary': 1.3e13, 'continuum': 1e14, 'line': 1e15},
      {'handler': 'move_packet_across_shell_boundary', 'distance': 1.3e13}),

     ({'boundary': 1.3e13, 'continuum': 1e14, 'line': 2.5e12},
      {'handler': 'montecarlo_line_scatter', 'distance': 2.5e12}),

     ({'boundary': 1.3e13, 'continuum': 1e11, 'line': 2.5e12},
      {'handler': 'montecarlo_thomson_scatter', 'distance': 1e11})]
)
def test_get_event_handler(packet, model, mt_state, distances, expected):
    d_cont_setter(distances['continuum'], model, packet)
    d_line_setter(distances['line'], model, packet)
    d_boundary_setter(distances['boundary'], model, packet)
    obtained_distance = c_double()

    cmontecarlo_methods.get_event_handler.restype = c_void_p
    obtained_handler = cmontecarlo_methods.get_event_handler(byref(packet), byref(model),
                                                             byref(obtained_distance),
                                                             byref(mt_state))

    expected_handler = getattr(cmontecarlo_methods, expected['handler'])
    expected_handler = cast(expected_handler, c_void_p).value

    assert_equal(obtained_handler, expected_handler)
    assert_allclose(obtained_distance.value, expected['distance'], rtol=1e-10)


"""
Simple Tests:
----------------
These test check very simple pices of code still work.
"""

@pytest.mark.parametrize(
    ['packet_params', 'line_idx', 'expected'],
    [({'energy': 0.0}, 0, 0),
     ({'energy': 1.0}, 1, 1),
     ({'energy': 0.5}, 2, 1.5)]
)
def test_increment_Edotlu_estimator(packet_params, line_idx, expected, packet, model):
    packet.energy = packet_params['energy']

    cmontecarlo_methods.increment_Edotlu_estimator(byref(packet), byref(model), c_int64(line_idx))

    assert_almost_equal(model.line_lists_Edotlu[line_idx], expected)

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



"""
Continuum Tests:
----------------
The tests written further (till next block comment is encountered) are for the
methods related to continuum interactions.
"""


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    't_electron', [2500., 15000.]
)
def test_sample_nu_free_free(t_electron, packet, model, mt_state_seeded, expected_ff_emissivity):
    model.t_electrons[packet.current_shell_id] = t_electron
    cmontecarlo_methods.sample_nu_free_free.restype = c_double

    nu_bins, expected_emissivity = expected_ff_emissivity(t_electron)

    nus = []
    for _ in xrange(int(1e5)):
        nu = cmontecarlo_methods.sample_nu_free_free(byref(packet), byref(model), byref(mt_state_seeded))
        nus.append(nu)

    obtained_emissivity, _ = np.histogram(nus, normed=True, bins=nu_bins)

    assert_equal(obtained_emissivity, expected_emissivity)


@pytest.mark.parametrize(
    ['z_random', 'current_shell_id', 'expected'],
    [(0.22443743797312765, 0, 0),
     (0.54510721066252377, 0, 2),
     (0.94183547596539363, 0, 3),
     (0.22443743797312765, 1, 2),
     (0.94183547596539363, 1, 3)]
)
@pytest.mark.continuumtest
def test_sample_cooling_processes(packet, z_random, current_shell_id, get_rkstate, expected):
    cooling_references = (c_int64 * 4)(0, 1, 2, 3)
    no_processes = c_int64(len(cooling_references))
    cooling_probabilities = (c_double * 8)(*([0.25] * 4 + [0.10] * 2 + [0.40] * 2))

    cmontecarlo_methods.sample_cooling_processes.restype = c_int64

    packet.current_shell_id = current_shell_id
    rkstate = get_rkstate(z_random)

    obtained_reference = cmontecarlo_methods.sample_cooling_processes(
        byref(packet), byref(rkstate), cooling_probabilities, cooling_references, no_processes)

    assert_equal(obtained_reference, expected)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['packet_params', 't_electrons', 'chi_ff_factor', 'expected'],
    [({'nu': 4.5e14, 'mu': 0.0, 'current_shell_id': 1}, 15000, 2.0, 1.6746639430359494e-44),
     ({'nu': 3.0e15, 'mu': 0.0, 'current_shell_id': 0}, 5000, 3.0, 1.1111111111107644e-46),
     ({'nu': 3.0e15, 'mu': 0.4, 'current_shell_id': 0}, 10000, 4.0, 1.5638286016098277e-46)]
)
def test_calculate_chi_ff(packet, model, packet_params, t_electrons, chi_ff_factor, expected):
    packet.mu = packet_params['mu']
    packet.nu = packet_params['nu']
    packet.current_shell_id = packet_params['current_shell_id']
    packet.r = 1.04e17

    model.t_electrons[packet_params['current_shell_id']] = t_electrons
    model.chi_ff_factor[packet_params['current_shell_id']] = chi_ff_factor

    cmontecarlo_methods.calculate_chi_ff(byref(packet), byref(model))
    obtained = packet.chi_ff

    assert_equal(obtained, expected)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['continuum_status', 'z_random', 'packet_params', 'expected'],
    [(CONTINUUM_OFF, 0.94183547596539363,
      {'chi_c': 1.0, 'chi_th': 0.4, 'chi_bf': 0.5},
      'montecarlo_thomson_scatter'),

     (CONTINUUM_ON, 0.22443743797312765,
      {'chi_c': 1.0, 'chi_th': 0.4, 'chi_bf': 0.5},
      'montecarlo_thomson_scatter'),

     (CONTINUUM_ON, 0.54510721066252377,
      {'chi_c': 1.0, 'chi_th': 0.4, 'chi_bf': 0.5},
      'montecarlo_bound_free_scatter'),

     (CONTINUUM_ON, 0.94183547596539363,
      {'chi_c': 1.0, 'chi_th': 0.4, 'chi_bf': 0.5},
      'montecarlo_free_free_scatter'),

     (CONTINUUM_ON, 0.22443743797312765,
      {'chi_c': 1e2, 'chi_th': 1e1, 'chi_bf': 2e1},
      'montecarlo_bound_free_scatter')]
)
def test_montecarlo_continuum_event_handler(continuum_status, expected, z_random,
                                            packet_params, packet, model, get_rkstate):
    packet.chi_cont = packet_params['chi_c']
    packet.chi_th = packet_params['chi_th']
    packet.chi_bf = packet_params['chi_bf']
    model.cont_status = continuum_status

    rkstate = get_rkstate(z_random)

    cmontecarlo_methods.montecarlo_continuum_event_handler.restype = c_void_p
    obtained = cmontecarlo_methods.montecarlo_continuum_event_handler(byref(packet),
                                                                      byref(model), byref(rkstate))
    expected = getattr(cmontecarlo_methods, expected)
    expected = cast(expected, c_void_p).value

    assert_equal(obtained, expected)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['nu', 'continuum_id', 'expected'],
    [(4.40e14, 1, 0.00),
     (3.25e14, 1, 0.75),
     (4.03e14, 0, 0.97),
     (4.10e14 + 1e-1, 0, 0.90),
     pytest.mark.xfail(reason="nu coincides with a supporting point")((4.1e14, 0, 0.90))]
)
def test_bf_cross_section(nu, continuum_id, model_w_edges, expected):
    cmontecarlo_methods.bf_cross_section.restype = c_double
    obtained = cmontecarlo_methods.bf_cross_section(byref(model_w_edges), continuum_id, c_double(nu))

    assert_almost_equal(obtained, expected)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['packet_params', 'expected'],
    [({'nu': 4.13e14, 'mu': 0.0, 'current_shell_id': 1},
      [3.2882087455641473] * 3),

     ({'nu': 3.27e14, 'mu': 0.0, 'current_shell_id': 0},
      [0.0, 1.3992114634681028, 5.702548202131454]),

     ({'nu': 3.27e14, 'mu': -0.4, 'current_shell_id': 0},
      [0.0, 1.2670858, 5.4446587])]
)
def test_calculate_chi_bf(packet_params, expected, packet, model_w_edges):
    model_w_edges.l_pop = (c_double * 6)(*range(1, 7))
    model_w_edges.l_pop_r = (c_double * 6)(*np.linspace(0.1, 0.6, 6))
    model_w_edges.t_electrons[packet_params['current_shell_id']] = 1e4

    packet.mu = packet_params['mu']
    packet.nu = packet_params['nu']
    packet.r = 1.04e17
    packet.current_shell_id = packet_params['current_shell_id']
    packet.chi_bf_tmp_partial = (c_double * model_w_edges.no_of_edges)()

    cmontecarlo_methods.calculate_chi_bf(byref(packet), byref(model_w_edges))

    obtained_chi_bf_tmp = np.ctypeslib.as_array(packet.chi_bf_tmp_partial, shape=(model_w_edges.no_of_edges,))
    expected_chi_bf_tmp = np.array(expected)
    expected_chi_bf = expected_chi_bf_tmp[-1]

    assert_almost_equal(obtained_chi_bf_tmp, expected_chi_bf_tmp)
    assert_almost_equal(packet.chi_bf, expected_chi_bf)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['comov_energy', 'distance', 'chi_ff', 'no_of_updates', 'expected'],
    [(1.3, 1.3e14, 3e-12, 1, 507.),
     (0.9, 0.7e15, 2e-12, 10, 1.260e4),
     (0.8, 0.8e13, 1.5e-12, 35, 336.)]
)
def test_increment_continuum_estimators_ff_heating_estimator(packet, model_w_edges, comov_energy, distance,
                                                             chi_ff, no_of_updates, expected):
    packet.chi_ff = chi_ff

    for _ in range(no_of_updates):
        cmontecarlo_methods.increment_continuum_estimators(byref(packet), byref(model_w_edges), c_double(distance),
                                                           c_double(0), c_double(comov_energy))
    obtained = model_w_edges.ff_heating_estimator[packet.current_shell_id]

    assert_almost_equal(obtained, expected)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['comov_nus', 'expected'],
    [([4.05e14, 4.17e14, 3.3e14, 3.2e14, 2.9e14], [2, 2, 3]),
     ([4.15e15, 3.25e14, 3.3e14, 2.85e14, 2.9e14], [0, 2, 4])]
)
def test_increment_continuum_estimators_photo_ion_estimator_statistics(packet, model_w_edges, comov_nus, expected):
    for comov_nu in comov_nus:
        cmontecarlo_methods.increment_continuum_estimators(byref(packet), byref(model_w_edges), c_double(1e13),
                                                           c_double(comov_nu), c_double(1.0))

    no_of_edges = model_w_edges.no_of_edges
    no_of_shells = model_w_edges.no_of_shells

    obtained = np.ctypeslib.as_array(model_w_edges.photo_ion_estimator_statistics,
                                     shape=(no_of_edges * no_of_shells,))
    obtained = np.reshape(obtained, newshape=(no_of_shells, no_of_edges), order='F')
    obtained = obtained[packet.current_shell_id]
    expected = np.array(expected)

    assert_array_equal(obtained, expected)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['comov_energy', 'distance', 'comov_nus', 'expected'],
    [(1.3, 1.3e14, [4.05e14, 2.65e14],
      {"photo_ion": [0.39641975308641975, 0., 0.],
       "stim_recomb": [0.056757061269242064, 0., 0.],
       "bf_heating": [1.9820987654321076e12, 0., 0.],
       "stim_recomb_cooling": [283784812699.75476, 0., 0.]}),

     (0.9, 0.7e15, [3.25e14, 2.85e14],
      {"photo_ion": [0., 1.4538461538461538, 7.315141700404858],
       "stim_recomb": [0., 0.3055802, 1.7292954],
       "bf_heating": [0., 36346153846153.82, 156760323886639.69],
       "stim_recomb_cooling": [0., 7639505724285.9746, 33907776077426.875]})]
)
def test_increment_continuum_estimators_bf_estimators(packet, model_w_edges, comov_energy,
                                                      distance, comov_nus, expected):
    for comov_nu in comov_nus:
        cmontecarlo_methods.increment_continuum_estimators(byref(packet), byref(model_w_edges), c_double(distance),
                                                           c_double(comov_nu), c_double(comov_energy))

    no_of_edges = model_w_edges.no_of_edges
    no_of_shells = model_w_edges.no_of_shells

    for estim_name, expected_value in expected.iteritems():
        obtained = np.ctypeslib.as_array(getattr(model_w_edges, estim_name + "_estimator"),
                                         shape=(no_of_edges * no_of_shells,))
        obtained = np.reshape(obtained, newshape=(no_of_shells, no_of_edges), order='F')
        obtained = obtained[packet.current_shell_id]

        assert_almost_equal(obtained, np.array(expected_value))


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['z_random', 'expected'],
    [(0.22443743797312765,  # next z_random = 0.545678896748
      {'process_type': BF_EMISSION, 'current_continuum_id': 2}),

     (0.58269284658380704,  # next z_random = 0.845102049919
      {'process_type': EXCITATION, 'macro_atom_activation_level': 3}),

     (0.78961460371187597,  # next z_random = 0.818455414618
      {'process_type': IONIZATION, 'macro_atom_activation_level': 3}),

     (0.3956855427967142,
      {'process_type': FF_EMISSION})]
)
def test_k_packet(packet, model_w_edges, get_rkstate, z_random, expected):
    rkstate = get_rkstate(z_random)

    model_w_edges.fb_cooling_prob[packet.current_shell_id] = 0.25
    model_w_edges.ff_cooling_prob[packet.current_shell_id] = 0.25
    model_w_edges.coll_exc_cooling_prob[packet.current_shell_id] = 0.25

    cooling_references = (c_int64 * 4)(0, 1, 2, 3)
    no_processes = c_int64(len(cooling_references))
    cooling_probabilities = (c_double * 8)(*([0.25] * 4 + [0.10] * 2 + [0.40] * 2))
    for name in ['fb', 'coll_exc', 'coll_ion']:
        setattr(model_w_edges, name + '_cooling_references', cooling_references)
        setattr(model_w_edges, name + '_cooling_prob_nd', no_processes)
        setattr(model_w_edges, name + '_cooling_prob_individual', cooling_probabilities)

    obtained_process_type = c_int(-10)
    expected_process_type = expected.pop('process_type')

    cmontecarlo_methods.k_packet(byref(packet), byref(model_w_edges), byref(obtained_process_type), byref(rkstate))

    assert_equal(obtained_process_type.value, expected_process_type)

    for reference_name, expected_reference_value in expected.iteritems():
        obtained_reference_value = getattr(packet, reference_name)
        assert_equal(obtained_reference_value, expected_reference_value)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    [({'nu_comov': 1.1e16, 'mu': 0.0, 'r': 1.4e14},
      {'next_line_id': 5, 'last_line': True, 'nu': 1.1e16}),

     ({'nu_comov': 1.3e16, 'mu': 0.3, 'r': 7.5e14},
      {'next_line_id': 0, 'last_line': False, 'nu': 1.30018766e+16}),

     ({'nu_comov': 1.24e16, 'mu': -0.3, 'r': 7.5e14},
      {'next_line_id': 2, 'last_line': False, 'nu': 1.23982106e+16})]
)
def test_continuum_emission(packet, model, mock_sample_nu, packet_params, expected_params, mt_state):
    packet.nu = packet_params['nu_comov']  # Is returned by mock function mock_sample_nu
    packet.mu = packet_params['mu']
    packet.r = packet_params['r']

    cmontecarlo_methods.continuum_emission(byref(packet), byref(model), byref(mt_state), mock_sample_nu)

    obtained_next_line_id = packet.next_line_id

    assert_equal(obtained_next_line_id, expected_params['next_line_id'])
    assert_equal(packet.last_line, expected_params['last_line'])
    assert_allclose(packet.nu, expected_params['nu'], rtol=1e-7)


"""
Not Yet Relevant Tests:
-----------------------
The tests written further (till next block comment is encountered) are for the
methods related to Continuum interactions. These are not required to be tested
on current master and can be skipped for now.
"""


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
