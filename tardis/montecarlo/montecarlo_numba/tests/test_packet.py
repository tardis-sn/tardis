import os
import pytest
import numpy as np
import pandas as pd
import tardis.montecarlo.formal_integral as formal_integral
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_configuration as mc
import tardis.montecarlo.montecarlo_numba.numba_interface as numba_interface
from tardis import constants as const
from tardis.montecarlo.montecarlo_numba.numba_interface import Estimators
from tardis.montecarlo.montecarlo_numba import macro_atom
C_SPEED_OF_LIGHT = const.c.to('cm/s').value

from numpy.testing import (
        assert_equal,
        assert_almost_equal,
        assert_array_equal,
        assert_allclose
        )

@pytest.fixture(scope="function")
def packet():
    return r_packet.RPacket(
        r = 7.5e14,
        nu = 0.4,
        mu = 0.3,
        energy = 0.9,
        seed = 1963,
        index = 0,
        is_close_line = 0
    )

@pytest.fixture(scope="function")
def model():
    return numba_interface.NumbaModel(
        r_inner = np.array([6.912e14, 8.64e14], dtype=np.float64),
        r_outer = np.array([8.64e14, 1.0368e15], dtype=np.float64),
        time_explosion = 5.2e7
    )


@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    [({'mu': 0.3, 'r': 7.5e14},
      {'d_boundary': 259376919351035.88}),

     ({'mu': -.3, 'r': 7.5e13},
      {'d_boundary': -664987228972291.5}),

     ({'mu': -.3, 'r': 7.5e14},
      {'d_boundary': 709376919351035.9})]
)
def test_calculate_distance_boundary(packet_params, expected_params, model):
    mu = packet_params['mu']
    r = packet_params['r']

    d_boundary = r_packet.calculate_distance_boundary(
        r, mu, model.r_inner[0], model.r_outer[0])

    #Accuracy to within 0.1cm
    assert_almost_equal(d_boundary[0], expected_params['d_boundary'], decimal=1)
#
#
# TODO: split this into two tests - one to assert errors and other for d_line
@pytest.mark.parametrize(
    ['packet_params', 'expected_params'],
    [({'nu_line': 0.1, 'next_line_id': 0, 'is_last_line': True},
      {'tardis_error': None, 'd_line': 1e+99}),

     ({'nu_line': 0.2, 'next_line_id': 1, 'is_last_line': False},
      {'tardis_error': None, 'd_line': 7.792353908000001e+17}),

     ({'nu_line': 0.5, 'next_line_id': 1, 'is_last_line': False},
      {'tardis_error': r_packet.MonteCarloException, 'd_line': 0.0}),

     ({'nu_line': 0.6, 'next_line_id': 0, 'is_last_line': False},
      {'tardis_error': r_packet.MonteCarloException, 'd_line': 0.0})]
)
def test_calculate_distance_line(packet_params, expected_params, packet, model):
    nu_line = packet_params['nu_line']
    is_last_line = packet_params['is_last_line']

    time_explosion = model.time_explosion

    doppler_factor = r_packet.get_doppler_factor(packet.r,
                                        packet.mu,
                                        time_explosion)
    comov_nu = packet.nu * doppler_factor

    d_line = 0
    obtained_tardis_error = None
    try:
        d_line = r_packet.calculate_distance_line(packet,
                                                    comov_nu,
                                                    is_last_line,
                                                    nu_line,
                                                    time_explosion)
    except r_packet.MonteCarloException:
        obtained_tardis_error = r_packet.MonteCarloException

    assert_almost_equal(d_line, expected_params['d_line'])
    assert obtained_tardis_error == expected_params['tardis_error']

def test_calculate_distance_electron():
    pass

def test_calculate_tau_electron():
    pass

@pytest.mark.parametrize(
    ['mu', 'r', 'inv_t_exp', 'expected'],
    [(0.3, 7.5e14, 1 / 5.2e7, 0.9998556693818854),
     (-.3, 0, 1 / 2.6e7, 1.0),
     (0, 1, 1 / 2.6e7, 1.0)]
)
def test_get_doppler_factor(mu, r, inv_t_exp, expected):
    # Set the params from test cases here
    # TODO: add relativity tests
    time_explosion = 1/inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = r_packet.get_doppler_factor(r, mu, time_explosion)


    # Perform required assertions
    assert_almost_equal(obtained, expected)

@pytest.mark.parametrize(
    ['mu', 'r', 'inv_t_exp', 'expected'],
    [(0.3, 7.5e14, 1 / 5.2e7, 1/0.9998556693818854),
     (-.3, 0, 1 / 2.6e7, 1.0),
     (0, 1, 1 / 2.6e7, 1.0)]
)
def test_get_inverse_doppler_factor(mu, r, inv_t_exp, expected):
    # Set the params from test cases here
    # TODO: add relativity tests
    time_explosion = 1/inv_t_exp

    # Perform any other setups just before this, they can be additional calls
    # to other methods or introduction of some temporary variables

    obtained = r_packet.get_inverse_doppler_factor(r, mu, time_explosion)

    # Perform required assertions
    assert_almost_equal(obtained, expected)


def test_get_random_mu():
    """
    Ensure that different calls results
    """
    output1 = r_packet.get_random_mu()
    output2 = r_packet.get_random_mu()
    assert output1 != output2

def test_update_line_estimators():
    pass

def test_trace_packet():
    pass

def test_move_r_packet():
    pass

def test_set_estimators():
    pass

def test_set_estimators_full_relativity():
    pass

def test_line_emission():
    pass

@pytest.mark.parametrize(
    ['current_shell_id', 'delta_shell', 'no_of_shells'],
    [(132, 11, 132),
     (132, 1, 133),
     (132, 2, 133)]
)
def test_move_packet_across_shell_boundary_emitted(packet, current_shell_id,
                                                   delta_shell,
                                                   no_of_shells):
    packet.current_shell_id = current_shell_id
    r_packet.move_packet_across_shell_boundary(packet, delta_shell,
                                      no_of_shells)
    assert packet.status == r_packet.PacketStatus.EMITTED

@pytest.mark.parametrize(
    ['current_shell_id', 'delta_shell', 'no_of_shells'],
    [(132, -133, 132),
     (132, -133, 133),
     (132, -1e9, 133)]
)
def test_move_packet_across_shell_boundary_reabsorbed(packet, current_shell_id,
                                                   delta_shell,
                                                   no_of_shells):
    packet.current_shell_id = current_shell_id
    r_packet.move_packet_across_shell_boundary(packet, delta_shell,
                                      no_of_shells)
    assert packet.status == r_packet.PacketStatus.REABSORBED


@pytest.mark.parametrize(
    ['current_shell_id', 'delta_shell', 'no_of_shells'],
    [(132, -1, 199),
     (132, 0, 132),
     (132, 20, 154)]
)
def test_move_packet_across_shell_boundary_increment(packet, current_shell_id,
                                                   delta_shell,
                                                   no_of_shells):
    packet.current_shell_id = current_shell_id
    r_packet.move_packet_across_shell_boundary(packet, delta_shell,
                                      no_of_shells)
    assert packet.current_shell_id == current_shell_id + delta_shell

#SAYS WE NEED TO FIX/REMOVE PACKET CALC ENERGY BY MOVING DOPPLER FACTOR TO FUNCTION
"""
@pytest.mark.parametrize(
    ['distance_trace', 'time_explosion'],
    [(0, 1),
     (1, 1),
     (1, 1e9)]
)
def test_packet_energy_limit_one(packet, distance_trace, time_explosion):
    initial_energy = packet.energy
    new_energy = r_packet.calc_packet_energy(packet, distance_trace, time_explosion)
    assert_almost_equal(new_energy, initial_energy)
"""
def test_test_for_close_line():
    pass