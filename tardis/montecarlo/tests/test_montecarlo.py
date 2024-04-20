import os

import numpy as np
import pandas as pd
import pytest

import tardis.montecarlo.montecarlo_numba.formal_integral as formal_integral
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_numba.utils as utils
import tardis.transport.r_packet_transport as r_packet_transport
from tardis import constants as const
from tardis.montecarlo.estimators.radfield_mc_estimators import (
    RadiationFieldMCEstimators,
)
from tardis.montecarlo.montecarlo_numba.numba_interface import RPacketTracker
from tardis.transport.frame_transformations import (
    angle_aberration_CMF_to_LF,
    angle_aberration_LF_to_CMF,
    get_doppler_factor,
)

C_SPEED_OF_LIGHT = const.c.to("cm/s").value

pytestmark = pytest.mark.skip(reason="Port from C to numba")


from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
)

from tardis import __path__ as path


@pytest.fixture(scope="module")
def continuum_compare_data_fname():
    fname = "continuum_compare_data.hdf"
    return os.path.join(path[0], "montecarlo", "tests", "data", fname)


@pytest.fixture(scope="module")
def continuum_compare_data(continuum_compare_data_fname, request):
    compare_data = pd.HDFStore(continuum_compare_data_fname, mode="r")

    def fin():
        compare_data.close()

    request.addfinalizer(fin)

    return compare_data


@pytest.fixture(scope="function")
def expected_ff_emissivity(continuum_compare_data):
    emissivities = continuum_compare_data["ff_emissivity"]

    def ff_emissivity(t_electron):
        emissivity = emissivities[t_electron]
        nu_bins = emissivity["nu_bins"].values
        emissivity_value = emissivity["emissivity"].dropna().values

        return nu_bins, emissivity_value

    return ff_emissivity


@pytest.fixture(scope="module")
def ion_edges():
    return [
        {
            "nu": [4.0e14, 4.1e14, 4.2e14, 4.3e14],
            "x_sect": [1.0, 0.9, 0.8, 0.7],
            "no_of_points": 4,
        },
        {
            "nu": [3.0e14, 3.1e14, 3.2e14, 3.3e14, 3.4e14],
            "x_sect": [1.0, 0.9, 0.8, 0.7, 0.6],
            "no_of_points": 5,
        },
        {
            "nu": [2.8e14, 3.0e14, 3.2e14, 3.4e14],
            "x_sect": [2.0, 1.8, 1.6, 1.4],
            "no_of_points": 4,
        },
    ]


"""
Important Tests:
----------------
The tests written further (till next block comment is encountered) have been
categorized as important tests, these tests correspond to methods which are
relatively old and stable code.
"""


@pytest.mark.parametrize(
    ["x", "x_insert", "imin", "imax", "expected_params"],
    [
        (
            [5.0, 4.0, 3.0, 1.0],
            2.0,
            0,
            3,
            {"result": 2},
        ),
        (
            [5.0, 4.0, 3.0, 2.0],
            0.0,
            0,
            3,
            {"result": 0},  # This one might need to check for a bounds error
        ),
    ],
)
def test_reverse_binary_search(x, x_insert, imin, imax, expected_params):
    # x = (c_double * (imax - imin + 1))(*x)
    obtained_result = 0

    obtained_result = formal_integral.reverse_binary_search(
        x, x_insert, imin, imax
    )

    assert obtained_result == expected_params["result"]


@pytest.mark.parametrize(
    ["nu", "nu_insert", "number_of_lines", "expected_params"],
    [
        (
            [0.5, 0.4, 0.3, 0.1],
            0.2,
            4,
            {"result": 3},
        ),
        (
            [0.5, 0.4, 0.3, 0.2],
            0.1,
            4,
            {"result": 4},
        ),
        (
            [0.4, 0.3, 0.2, 0.1],
            0.5,
            4,
            {"result": 0},
        ),
    ],
)
def test_line_search(nu, nu_insert, number_of_lines, expected_params):
    # nu = (c_double * number_of_lines)(*nu)
    obtained_result = 0

    obtained_result = formal_integral.line_search(
        nu, nu_insert, number_of_lines, obtained_result
    )

    assert obtained_result == expected_params["result"]


@pytest.mark.parametrize(
    ["x", "x_insert", "imin", "imax", "expected_params"],
    [
        (
            [2.0, 4.0, 6.0, 7.0],
            5.0,
            0,
            3,
            {"result": 2},
        ),
        (
            [2.0, 3.0, 5.0, 7.0],
            8.0,
            0,
            3,
            {"result": 0},  # this one too
        ),
        (
            [2.0, 4.0, 6.0, 7.0],
            4.0,
            0,
            3,
            {"result": 0},
        ),
    ],
)
def test_binary_search(x, x_insert, imin, imax, expected_params):
    obtained_result = 0

    obtained_result = formal_integral.binary_search(x, x_insert, imin, imax)

    assert obtained_result == expected_params["result"]


def test_get_random_mu_different_output():
    """
    Ensure that different calls results
    """
    output1 = r_packet.get_random_mu()
    output2 = r_packet.get_random_mu()
    assert output1 != output2


def test_get_random_mu_different_output():
    """
    Ensure that different calls results
    """
    output1 = r_packet.get_random_mu()
    output2 = r_packet.get_random_mu()
    assert output1 != output2


@pytest.mark.parametrize(
    ["mu", "r", "time_explosion"], [(1, C_SPEED_OF_LIGHT, 1)]
)
def test_angle_ab_LF_to_CMF_diverge(mu, r, time_explosion):
    """
    This equation should diverage as costheta --> 1 and beta --> 1
    """
    nu = 0.4
    energy = 0.9
    packet = r_packet.RPacket(r, mu, nu, energy)
    with pytest.raises(ZeroDivisionError):
        obtained = r_packet.angle_aberration_LF_to_CMF(
            packet, time_explosion, mu
        )


@pytest.mark.parametrize(
    ["mu", "r", "time_explosion"], [(0.3, 1e7, 1e10), (-0.3, 1e7, 1e11)]
)
def test_both_angle_aberrations(mu, r, time_explosion):
    """
    The angle aberration functions should be the functional inverse of one
    another.
    """
    nu = 0.4
    energy = 0.9
    packet = r_packet.RPacket(r, mu, nu, energy)
    packet.r = r
    obtained_mu = angle_aberration_LF_to_CMF(packet, time_explosion, mu)
    inverse_obtained_mu = angle_aberration_CMF_to_LF(
        packet, time_explosion, obtained_mu
    )
    assert_almost_equal(inverse_obtained_mu, mu)


@pytest.mark.parametrize(
    ["mu", "r", "time_explosion"], [(0.3, 7.5e14, 5.2e5), (-0.3, 7.5e14, 5.2e5)]
)
def test_both_angle_aberrations_inverse(mu, r, time_explosion):
    """
    The angle aberration functions should be the functional inverse of one
    another.
    """
    nu = 0.4
    energy = 0.9
    packet = r_packet.RPacket(r, mu, nu, energy)
    packet.r = r
    obtained_mu = angle_aberration_CMF_to_LF(packet, time_explosion, mu)
    inverse_obtained_mu = angle_aberration_LF_to_CMF(
        packet, time_explosion, obtained_mu
    )
    assert_almost_equal(inverse_obtained_mu, mu)


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, 11, 132), (132, 1, 133), (132, 2, 133)],
)
def test_move_packet_across_shell_boundary_emitted(
    current_shell_id, delta_shell, no_of_shells
):
    r = 7.5e14
    mu = 0.3
    nu = 0.4
    energy = 0.9
    packet = r_packet.RPacket(r, mu, nu, energy)
    packet.current_shell_id = current_shell_id
    r_packet_transport.move_packet_across_shell_boundary(
        packet, delta_shell, no_of_shells
    )
    assert packet.status == r_packet.PacketStatus.EMITTED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -133, 132), (132, -133, 133), (132, -1e9, 133)],
)
def test_move_packet_across_shell_boundary_reabsorbed(
    current_shell_id, delta_shell, no_of_shells
):
    r = 7.5e14
    mu = 0.3
    nu = 0.4
    energy = 0.9
    packet = r_packet.RPacket(r, mu, nu, energy)
    packet.current_shell_id = current_shell_id
    r_packet_transport.move_packet_across_shell_boundary(
        packet, delta_shell, no_of_shells
    )
    assert packet.status == r_packet.PacketStatus.REABSORBED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -1, 199), (132, 0, 132), (132, 20, 154)],
)
def test_move_packet_across_shell_boundary_increment(
    current_shell_id, delta_shell, no_of_shells
):
    r = 7.5e14
    mu = 0.3
    nu = 0.4
    energy = 0.9
    packet = r_packet.RPacket(r, mu, nu, energy)
    packet.current_shell_id = current_shell_id
    r_packet_transport.move_packet_across_shell_boundary(
        packet, delta_shell, no_of_shells
    )
    assert packet.current_shell_id == current_shell_id + delta_shell


@pytest.mark.parametrize(
    ["distance_trace", "time_explosion", "mu", "r"],
    [(0, 1, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1)],
)
def test_packet_energy_limit_one(distance_trace, time_explosion, mu, r):
    initial_energy = 0.9
    nu = 0.4
    packet = r_packet.RPacket(r, mu, nu, initial_energy)
    new_energy = r_packet.calc_packet_energy(
        packet, distance_trace, time_explosion
    )
    assert new_energy == initial_energy


@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        ({"mu": 0.3, "r": 7.5e14}, {"d_boundary": 259376919351035.88}),
        ({"mu": -0.3, "r": 7.5e13}, {"d_boundary": -664987228972291.5}),
        ({"mu": -0.3, "r": 7.5e14}, {"d_boundary": 709376919351035.9}),
    ],
)
def test_compute_distance2boundary(packet_params, expected_params):
    mu = packet_params["mu"]
    r = packet_params["r"]
    r_inner = np.array([6.912e14, 8.64e14], dtype=np.float64)
    r_outer = np.array([8.64e14, 1.0368e15], dtype=np.float64)

    d_boundary = r_packet.calculate_distance_boundary(
        r, mu, r_inner[0], r_outer[0]
    )

    assert_almost_equal(d_boundary[0], expected_params["d_boundary"])


#
#
# TODO: split this into two tests - one to assert errors and other for d_line
@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        (
            {"nu_line": 0.1, "next_line_id": 0, "last_line": 1},
            {"tardis_error": None, "d_line": 1e99},
        ),
        (
            {"nu_line": 0.2, "next_line_id": 1, "last_line": 0},
            {"tardis_error": None, "d_line": 7.792353908000001e17},
        ),
        (
            {"nu_line": 0.5, "next_line_id": 1, "last_line": 0},
            {"tardis_error": utils.MonteCarloException, "d_line": 0.0},
        ),
        (
            {"nu_line": 0.6, "next_line_id": 0, "last_line": 0},
            {"tardis_error": utils.MonteCarloException, "d_line": 0.0},
        ),
    ],
)
def test_compute_distance2line(packet_params, expected_params):
    r = 7.5e14
    mu = 0.3
    nu = 0.4
    energy = 0.9
    packet = r_packet.RPacket(r, mu, nu, energy)
    nu_line = packet_params["nu_line"]
    # packet.next_line_id = packet_params['next_line_id']
    # packet.last_line = packet_params['last_line']

    time_explosion = 5.2e7

    doppler_factor = get_doppler_factor(
        packet.r, packet.mu, time_explosion, False
    )
    comov_nu = packet.nu * doppler_factor

    d_line = 0
    obtained_tardis_error = None
    try:
        d_line = r_packet.calculate_distance_line(
            packet, comov_nu, nu_line, time_explosion
        )
    except utils.MonteCarloException:
        obtained_tardis_error = utils.MonteCarloException

    assert_almost_equal(d_line, expected_params["d_line"])
    assert obtained_tardis_error == expected_params["tardis_error"]


# @pytest.mark.parametrize(
#     ['packet_params', 'expected_params'],
#     [({'virtual_packet': 0},
#      {'chi_cont': 6.652486e-16, 'd_cont': 4.359272608766106e+28}),
#
#      ({'virtual_packet': 1},
#       {'chi_cont': 6.652486e-16, 'd_cont': 1e+99})]
# )
# def test_compute_distance2continuum(clib, packet_params, expected_params, packet, model):
#     packet.virtual_packet = packet_params['virtual_packet']
#
#     clib.compute_distance2continuum(byref(packet), byref(model))
#
#     assert_almost_equal(packet.chi_cont, expected_params['chi_cont'])
#     assert_almost_equal(packet.d_cont, expected_params['d_cont'])


@pytest.mark.parametrize("full_relativity", [1, 0])
@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        (
            {"nu": 0.4, "mu": 0.3, "energy": 0.9, "r": 7.5e14},
            {
                "mu": 0.3120599529139568,
                "r": 753060422542573.9,
                "j": 8998701024436.969,
                "nubar": 3598960894542.354,
            },
        ),
        (
            {"nu": 0.6, "mu": -0.5, "energy": 0.5, "r": 8.1e14},
            {
                "mu": -0.4906548373534084,
                "r": 805046582503149.2,
                "j": 5001298975563.031,
                "nubar": 3001558973156.1387,
            },
        ),
    ],
)
def test_move_packet(packet_params, expected_params, full_relativity):
    distance = 1e13
    r, mu, nu, energy = 7.5e14, 0.3, 0.4, 0.9
    time_explosion = 5.2e7
    packet = r_packet.RPacket(r, mu, nu, energy)
    packet.nu = packet_params["nu"]
    packet.mu = packet_params["mu"]
    packet.energy = packet_params["energy"]
    packet.r = packet_params["r"]
    # model.full_relativity = full_relativity

    doppler_factor = get_doppler_factor(
        packet.r, packet.mu, time_explosion, full_relativity
    )
    numba_estimator = RadiationFieldMCEstimators(
        packet_params["j"], packet_params["nu_bar"], 0, 0
    )
    r_packet_transport.move_r_packet(
        packet, distance, time_explosion, numba_estimator, full_relativity
    )

    assert_almost_equal(packet.mu, expected_params["mu"])
    assert_almost_equal(packet.r, expected_params["r"])

    expected_j = expected_params["j"]
    expected_nubar = expected_params["nubar"]
    if full_relativity:
        expected_j *= doppler_factor
        expected_nubar *= doppler_factor

    mc.ENABLE_FULL_RELATIVITY = False

    assert_allclose(
        numba_estimator.j_estimator[packet.current_shell_id],
        expected_j,
        rtol=5e-7,
    )
    assert_allclose(
        numba_estimator.nu_bar_estimator[packet.current_shell_id],
        expected_nubar,
        rtol=5e-7,
    )


# @pytest.mark.continuumtest
# @pytest.mark.parametrize(
#     ['packet_params', 'j_blue_idx', 'expected'],
#     [({'nu': 0.30, 'energy': 0.30}, 0, 1.0),
#      ({'nu': 0.20, 'energy': 1.e5}, 0, 5e5),
#      ({'nu': 2e15, 'energy': 0.50}, 1, 2.5e-16),
#      ({'nu': 0.40, 'energy': 1e-7}, 1, 2.5e-7)],
# )
# def test_increment_j_blue_estimator_full_relativity(packet_params,
#                                                     j_blue_idx, expected,
#                                                     packet, model):
#     packet.nu = packet_params['nu']
#     packet.energy = packet_params['energy']
#     model.full_relativity = True
#
#     r_packet.increment_j_blue_estimator(byref(packet), byref(model),
#                                     c_double(packet.d_line),
#                                     c_int64(j_blue_idx))
#
#     assert_almost_equal(model.line_lists_j_blues[j_blue_idx], expected)
#
#
# @pytest.mark.parametrize(
#     ['packet_params', 'cur_line_id', 'expected'],
#     [({'nu': 0.1, 'mu': 0.3, 'r': 7.5e14}, 0, 8.998643292289723),
#      ({'nu': 0.2, 'mu': -.3, 'r': 7.7e14}, 0, 4.499971133976377),
#      ({'nu': 0.5, 'mu': 0.5, 'r': 7.9e14}, 1, 0.719988453650551),
#      ({'nu': 0.6, 'mu': -.5, 'r': 8.1e14}, 1, 0.499990378058792)]
# )
# def test_increment_j_blue_estimator(packet_params, cur_line_id, expected, packet):
#
#     numba_interface
#     packet = r_packet.RPacket(packet_params['r'],
#                               packet_params['mu'],
#                               packet_params['nu'],
#                               packet.energy)
#
#     r_packet.compute_distance2line(byref(packet), byref(model))
#     r_packet.move_r_packet(packet,
#                            distance,
#                            model.time_explosion,
#                            numba_estimator)
#     r_packet.move_packet(byref(packet), byref(model), c_double(1.e13))
#     r_packet.update_line_estimators(estimators, r_packet,
#                                     cur_line_id,
#                                     model.distance_trace,
#                                     model.time_explosion)
#
#     assert_almost_equal(model.line_lists_j_blues[j_blue_idx], expected)


"""
Continuum Tests:
----------------
The tests written further (till next block comment is encountered) are for the
methods related to continuum interactions.
"""


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp", "full_relativity"],
    [
        (0.8, 7.5e14, 1 / 5.2e5, 1),
        (-0.7, 7.5e14, 1 / 5.2e5, 1),
        (0.3, 7.5e14, 1 / 2.2e5, 1),
        (0.0, 7.5e14, 1 / 2.2e5, 1),
        (-0.7, 7.5e14, 1 / 5.2e5, 0),
    ],
)
def test_frame_transformations(mu, r, inv_t_exp, full_relativity):
    packet = r_packet.RPacket(r=r, mu=mu, energy=0.9, nu=0.4)
    mc.ENABLE_FULL_RELATIVITY = bool(full_relativity)
    mc.ENABLE_FULL_RELATIVITY = full_relativity

    inverse_doppler_factor = r_packet.get_inverse_doppler_factor(
        r, mu, 1 / inv_t_exp
    )
    r_packet.angle_aberration_CMF_to_LF(packet, 1 / inv_t_exp, packet.mu)

    doppler_factor = get_doppler_factor(r, mu, 1 / inv_t_exp)
    mc.ENABLE_FULL_RELATIVITY = False

    assert_almost_equal(doppler_factor * inverse_doppler_factor, 1.0)


@pytest.mark.continuumtest
@pytest.mark.parametrize(
    ["mu", "r", "inv_t_exp"],
    [
        (0.8, 7.5e14, 1 / 5.2e5),
        (-0.7, 7.5e14, 1 / 5.2e5),
        (0.3, 7.5e14, 1 / 2.2e5),
        (0.0, 7.5e14, 1 / 2.2e5),
        (-0.7, 7.5e14, 1 / 5.2e5),
    ],
)
def test_angle_transformation_invariance(mu, r, inv_t_exp):
    packet = r_packet.RPacket(r, mu, 0.4, 0.9)
    mc.ENABLE_FULL_RELATIVITY = True

    mu1 = angle_aberration_CMF_to_LF(packet, 1 / inv_t_exp, mu)
    mu_obtained = angle_aberration_LF_to_CMF(packet, 1 / inv_t_exp, mu1)

    mc.ENABLE_FULL_RELATIVITY = False
    assert_almost_equal(mu_obtained, mu)


@pytest.mark.continuumtest
@pytest.mark.parametrize("full_relativity", [1, 0])
@pytest.mark.parametrize(
    ["mu", "r", "t_exp", "nu", "nu_line"],
    [
        (0.8, 7.5e14, 5.2e5, 1.0e15, 9.4e14),
        (0.0, 6.3e14, 2.2e5, 6.0e12, 5.8e12),
        (1.0, 9.0e14, 2.2e5, 4.0e8, 3.4e8),
        (0.9, 9.0e14, 0.5e5, 1.0e15, 4.5e14),
        (-0.7, 7.5e14, 5.2e5, 1.0e15, 9.8e14),
        (-1.0, 6.3e14, 2.2e5, 6.0e12, 6.55e12),
    ],
)
def test_compute_distance2line_relativistic(
    mu, r, t_exp, nu, nu_line, full_relativity
):
    packet = r_packet.RPacket(r=r, nu=nu, mu=mu, energy=0.9)
    # packet.nu_line = nu_line
    numba_estimator = RadiationFieldMCEstimators(
        transport.j_estimator,
        transport.nu_bar_estimator,
        transport.j_blue_estimator,
        transport.Edotlu_estimator,
    )
    mc.ENABLE_FULL_RELATIVITY = bool(full_relativity)

    doppler_factor = get_doppler_factor(r, mu, t_exp)
    comov_nu = packet.nu * doppler_factor
    distance = r_packet.calculate_distance_line(
        packet, comov_nu, nu_line, t_exp
    )
    r_packet_transport.move_r_packet(
        packet, distance, t_exp, numba_estimator, bool(full_relativity)
    )

    doppler_factor = get_doppler_factor(r, mu, t_exp)
    comov_nu = packet.nu * doppler_factor
    mc.ENABLE_FULL_RELATIVITY = False

    assert_allclose(comov_nu, nu_line, rtol=1e-14)


"""
Tests for Tracking RPacket Properties
"""


@pytest.mark.parametrize("seed", [2274437677])
@pytest.mark.parametrize(
    ["index", "r", "nu", "mu", "energy"],
    [
        (0, 0, 0, 0, 0),
        (10, 1.43e15, 6.57e14, 0.92701038, 0.00010347),
        (20, 2.04e15, 9.59e14, 0.96005622, 0.00010347),
        (30, 1.12e15, 2.17e14, 0.96787487, 0.00010347),
        (40, 2.15e15, 9.17e14, 0.96428633, 0.00010347),
        (100000, 2.23e15, 8.56e14, 0.96946856, 0.00010347),
    ],
)
def test_rpacket_tracking(index, seed, r, nu, mu, energy):
    # Setup Montecarlo_Configuration.INITIAL_TRACKING_ARRAY_LENGTH
    mc.INITIAL_TRACKING_ARRAY_LENGTH = 10

    tracked_rpacket_properties = RPacketTracker()
    test_rpacket = r_packet.RPacket(
        index=index,
        seed=seed,
        r=r,
        nu=nu,
        mu=mu,
        energy=energy,
    )

    # TearDown Montecarlo_Configuration.INITIAL_TRACKING_ARRAY_LENGTH
    mc.INITIAL_TRACKING_ARRAY_LENGTH = None

    tracked_rpacket_properties.track(test_rpacket)
    tracked_rpacket_properties.finalize_array()

    assert test_rpacket.index == tracked_rpacket_properties.index
    assert test_rpacket.seed == tracked_rpacket_properties.seed
    assert test_rpacket.r == tracked_rpacket_properties.r
    assert test_rpacket.nu == tracked_rpacket_properties.nu
    assert test_rpacket.mu == tracked_rpacket_properties.mu
    assert test_rpacket.energy == tracked_rpacket_properties.energy
    assert tracked_rpacket_properties.num_interactions == 1
