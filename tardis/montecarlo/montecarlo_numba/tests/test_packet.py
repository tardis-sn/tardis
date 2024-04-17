import numpy as np
import pytest

import tardis.montecarlo.montecarlo_configuration as numba_config
import tardis.montecarlo.estimators.radfield_mc_estimators
import tardis.montecarlo.montecarlo_numba.numba_interface as numba_interface
import tardis.montecarlo.montecarlo_numba.opacities as opacities
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_numba.utils as utils
import tardis.transport.frame_transformations as frame_transformations
import tardis.transport.geometry.calculate_distances as calculate_distances
import tardis.transport.r_packet_transport as r_packet_transport
from tardis import constants as const
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.montecarlo.estimators.radfield_estimator_calcs import (
    update_line_estimators,
)

C_SPEED_OF_LIGHT = const.c.to("cm/s").value
SIGMA_THOMSON = const.sigma_T.to("cm^2").value

from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
)


@pytest.fixture(scope="function")
def geometry():
    return NumbaRadial1DGeometry(
        r_inner=np.array([6.912e14, 8.64e14], dtype=np.float64),
        r_outer=np.array([8.64e14, 1.0368e15], dtype=np.float64),
        v_inner=np.array([-1, -1], dtype=np.float64),
        v_outer=np.array([-1, -1], dtype=np.float64),
    )


@pytest.fixture(scope="function")
def model():
    return numba_interface.NumbaModel(
        time_explosion=5.2e7,
    )


@pytest.fixture(scope="function")
def estimators():
    return tardis.montecarlo.estimators.radfield_mc_estimators.RadiationFieldMCEstimators(
        j_estimator=np.array([0.0, 0.0], dtype=np.float64),
        nu_bar_estimator=np.array([0.0, 0.0], dtype=np.float64),
        j_blue_estimator=np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64
        ),
        Edotlu_estimator=np.array(
            [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], dtype=np.float64
        ),
        photo_ion_estimator=np.empty((0, 0), dtype=np.float64),
        stim_recomb_estimator=np.empty((0, 0), dtype=np.float64),
        bf_heating_estimator=np.empty((0, 0), dtype=np.float64),
        stim_recomb_cooling_estimator=np.empty((0, 0), dtype=np.float64),
        photo_ion_estimator_statistics=np.empty((0, 0), dtype=np.int64),
    )


@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        ({"mu": 0.3, "r": 7.5e14}, {"d_boundary": 259376919351035.88}),
        ({"mu": -0.3, "r": 7.5e13}, {"d_boundary": -664987228972291.5}),
        ({"mu": -0.3, "r": 7.5e14}, {"d_boundary": 709376919351035.9}),
    ],
)
def test_calculate_distance_boundary(packet_params, expected_params, geometry):
    mu = packet_params["mu"]
    r = packet_params["r"]

    d_boundary = calculate_distances.calculate_distance_boundary(
        r, mu, geometry.r_inner[0], geometry.r_outer[0]
    )

    # Accuracy to within 0.1cm
    assert_almost_equal(d_boundary[0], expected_params["d_boundary"], decimal=1)


#
#
# TODO: split this into two tests - one to assert errors and other for d_line
@pytest.mark.parametrize(
    ["packet_params", "expected_params"],
    [
        (
            {"nu_line": 0.1, "next_line_id": 0, "is_last_line": True},
            {"tardis_error": None, "d_line": 1e99},
        ),
        (
            {"nu_line": 0.2, "next_line_id": 1, "is_last_line": False},
            {"tardis_error": None, "d_line": 7.792353908000001e17},
        ),
        (
            {"nu_line": 0.5, "next_line_id": 1, "is_last_line": False},
            {"tardis_error": utils.MonteCarloException, "d_line": 0.0},
        ),
        (
            {"nu_line": 0.6, "next_line_id": 0, "is_last_line": False},
            {"tardis_error": utils.MonteCarloException, "d_line": 0.0},
        ),
    ],
)
def test_calculate_distance_line(
    packet_params, expected_params, static_packet, model
):
    nu_line = packet_params["nu_line"]
    is_last_line = packet_params["is_last_line"]

    time_explosion = model.time_explosion

    doppler_factor = frame_transformations.get_doppler_factor(
        static_packet.r, static_packet.mu, time_explosion, False
    )
    comov_nu = static_packet.nu * doppler_factor

    d_line = 0
    obtained_tardis_error = None
    try:
        d_line = calculate_distances.calculate_distance_line(
            static_packet,
            comov_nu,
            is_last_line,
            nu_line,
            time_explosion,
            enable_full_relativity=False,
        )
    except utils.MonteCarloException:
        obtained_tardis_error = utils.MonteCarloException

    assert_almost_equal(d_line, expected_params["d_line"])
    assert obtained_tardis_error == expected_params["tardis_error"]


@pytest.mark.parametrize(
    ["electron_density", "tau_event"], [(1e-5, 1.0), (1e10, 1e10)]
)
def test_calculate_distance_electron(electron_density, tau_event):
    actual = calculate_distances.calculate_distance_electron(
        electron_density, tau_event
    )
    expected = tau_event / (electron_density * SIGMA_THOMSON)

    assert_almost_equal(actual, expected)


@pytest.mark.parametrize(
    ["electron_density", "distance"],
    [(1e-5, 1.0), (1e10, 1e10), (-1, 0), (-1e10, -1e10)],
)
def test_calculate_tau_electron(electron_density, distance):
    actual = opacities.calculate_tau_electron(electron_density, distance)
    expected = electron_density * SIGMA_THOMSON * distance

    assert_almost_equal(actual, expected)


def test_get_random_mu(set_seed_fixture):
    """
    Ensure that different calls results
    """
    set_seed_fixture(1963)

    output1 = utils.get_random_mu()
    assert output1 == 0.9136407866175174


@pytest.mark.parametrize(
    [
        "cur_line_id",
        "distance_trace",
        "time_explosion",
        "expected_j_blue",
        "expected_Edotlu",
    ],
    [
        (
            0,
            1e12,
            5.2e7,
            [[2.249673812803061, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[2.249673812803061 * 0.4, 0.0, 1.0], [0.0, 0.0, 1.0]],
        ),
        (
            0,
            0,
            5.2e7,
            [[2.249675256109242, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [
                [2.249675256109242 * 0.4, 0.0, 1.0],
                [0.0, 0.0, 1.0],
            ],
        ),
        (
            1,
            1e5,
            1e10,
            [[0.0, 0.0, 0.0], [2.249998311331767, 0.0, 0.0]],
            [[0.0, 0.0, 1.0], [2.249998311331767 * 0.4, 0.0, 1.0]],
        ),
    ],
)
def test_update_line_estimators(
    estimators,
    static_packet,
    cur_line_id,
    distance_trace,
    time_explosion,
    expected_j_blue,
    expected_Edotlu,
):
    update_line_estimators(
        estimators,
        static_packet,
        cur_line_id,
        distance_trace,
        time_explosion,
        enable_full_relativity=False,
    )

    assert_allclose(estimators.j_blue_estimator, expected_j_blue)
    assert_allclose(estimators.Edotlu_estimator, expected_Edotlu)


# TODO set RNG consistently
# TODO: update this test to use the correct trace_packet
@pytest.mark.xfail(reason="Need to fix estimator differences across runs")
def test_trace_packet(
    packet,
    verysimple_numba_model,
    verysimple_opacity_state,
    verysimple_estimators,
    set_seed_fixture,
):
    set_seed_fixture(1963)
    packet.initialize_line_id(verysimple_opacity_state, verysimple_numba_model)
    distance, interaction_type, delta_shell = r_packet_transport.trace_packet(
        packet,
        verysimple_numba_model,
        verysimple_opacity_state,
        verysimple_estimators,
        continuum_processes_enabled=False,
        enable_full_relativity=False,
        disable_line_scattering=False,
    )

    assert delta_shell == 1
    assert interaction_type == 3
    assert_almost_equal(distance, 22978745222176.88)


@pytest.mark.xfail(reason="bug in full relativity")
@pytest.mark.parametrize("ENABLE_FULL_RELATIVITY", [True, False])
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
def test_move_r_packet(
    packet_params,
    expected_params,
    packet,
    model,
    estimators,
    ENABLE_FULL_RELATIVITY,
):
    distance = 1.0e13
    packet.nu = packet_params["nu"]
    packet.mu = packet_params["mu"]
    packet.energy = packet_params["energy"]
    packet.r = packet_params["r"]

    numba_config.ENABLE_FULL_RELATIVITY = ENABLE_FULL_RELATIVITY
    r_packet_transport.move_r_packet.recompile()  # This must be done as move_r_packet was jitted with ENABLE_FULL_RELATIVITY
    doppler_factor = frame_transformations.get_doppler_factor(
        packet.r, packet.mu, model.time_explosion, ENABLE_FULL_RELATIVITY
    )

    r_packet_transport.move_r_packet(
        packet,
        distance,
        model.time_explosion,
        estimators,
        ENABLE_FULL_RELATIVITY,
    )

    assert_almost_equal(packet.mu, expected_params["mu"])
    assert_almost_equal(packet.r, expected_params["r"])

    expected_j = expected_params["j"]
    expected_nubar = expected_params["nubar"]

    if ENABLE_FULL_RELATIVITY:
        expected_j *= doppler_factor
        expected_nubar *= doppler_factor

    numba_config.ENABLE_FULL_RELATIVITY = False
    assert_allclose(
        estimators.j_estimator[packet.current_shell_id], expected_j, rtol=5e-7
    )
    assert_allclose(
        estimators.nu_bar_estimator[packet.current_shell_id],
        expected_nubar,
        rtol=5e-7,
    )


@pytest.mark.xfail(reason="To be implemented")
def test_set_estimators():
    pass


@pytest.mark.xfail(reason="To be implemented")
def test_set_estimators_full_relativity():
    pass


@pytest.mark.xfail(reason="To be implemented")
def test_line_emission():
    pass


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, 11, 132), (132, 1, 133), (132, 2, 133)],
)
def test_move_packet_across_shell_boundary_emitted(
    packet, current_shell_id, delta_shell, no_of_shells
):
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
    packet, current_shell_id, delta_shell, no_of_shells
):
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
    packet, current_shell_id, delta_shell, no_of_shells
):
    packet.current_shell_id = current_shell_id
    r_packet_transport.move_packet_across_shell_boundary(
        packet, delta_shell, no_of_shells
    )
    assert packet.current_shell_id == current_shell_id + delta_shell
