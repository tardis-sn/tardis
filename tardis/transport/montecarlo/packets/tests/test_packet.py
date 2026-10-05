import numpy as np
import pytest
from astropy import units as u
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
)

import tardis.opacities.opacities as opacities
import tardis.transport.frame_transformations as frame_transformations
import tardis.transport.geometry.calculate_distances as calculate_distances
import tardis.transport.montecarlo.configuration.montecarlo_globals as montecarlo_globals
import tardis.transport.montecarlo.modes.homologous_rad_packet_transport as r_packet_transport
import tardis.transport.montecarlo.packets.radiative_packet as radiative_packet
import tardis.transport.montecarlo.utils as utils
from tardis import constants as const
from tardis.model.geometry.radial1d_homologous import HomologousRadial1DGeometry
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_line,
)
from tardis.transport.montecarlo.packets.movement import (
    move_packet_across_shell_boundary,
    move_r_packet,
)
from tardis.transport.montecarlo.packets.radiative_packet import InteractionType

C_SPEED_OF_LIGHT = const.c.to("cm/s").value
SIGMA_THOMSON = const.sigma_T.to("cm^2").value

# Shell used by the boundary-distance tests; the same radii as shell 0 of the
# ``geometry`` fixture below [cm].
SHELL_R_INNER = 6.912e14
SHELL_R_OUTER = 8.64e14
# Absolute tolerance for distances [cm]. One ulp at r ~ 1e15 cm is 0.125 cm and
# the boundary roots combine a handful of operations on squared radii, so a few
# ulp of rounding is expected; 1 cm is a small multiple of that and is still
# 14 orders of magnitude below the shell width.
DISTANCE_ATOL_CM = 1.0
# Relative tolerance on the radius a packet lands at after travelling the
# returned distance. Double-precision eps is 2.2e-16; 1e-12 leaves headroom for
# fastmath reassociation and the sqrt/square round trip, while still rejecting
# any algebraic error in the root (which shifts the landing point by O(1)).
LANDING_RTOL = 1e-12


def propagate_along_ray(
    r: float, mu: float, distance: float
) -> tuple[float, float]:
    """
    Advance a straight ray in spherical geometry by the law of cosines.

    This is the forward map that the boundary and line distance solvers invert,
    so it gives an independent check of their roots.

    Parameters
    ----------
    r : float
        Starting radius [cm].
    mu : float
        Starting direction cosine relative to the radial direction.
    distance : float
        Path length travelled [cm].

    Returns
    -------
    r_new : float
        Radius after travelling ``distance`` [cm].
    mu_new : float
        Direction cosine at the new position (radial projection of the
        unchanged direction vector).
    """
    r_new = np.sqrt(r * r + distance * distance + 2.0 * r * distance * mu)
    mu_new = (r * mu + distance) / r_new
    return r_new, mu_new


@pytest.fixture(scope="function")
def geometry():
    time_explosion = 5.2e7
    r_inner = np.array([6.912e14, 8.64e14], dtype=np.float64)
    r_outer = np.array([8.64e14, 1.0368e15], dtype=np.float64)
    time_explosion_quantity = time_explosion * u.s
    return HomologousRadial1DGeometry(
        v_inner=r_inner * u.cm / time_explosion_quantity,
        v_outer=r_outer * u.cm / time_explosion_quantity,
        v_inner_boundary=None,
        v_outer_boundary=None,
        time_explosion=time_explosion_quantity,
    ).to_numba()


@pytest.fixture(scope="function")
def time_explosion():
    return 5.2e7


@pytest.fixture(scope="function")
def estimators():
    from tardis.transport.montecarlo.estimators import EstimatorsLine

    return EstimatorsLine(
        mean_intensity_blueward=np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64
        ),
        energy_deposition_line_rate=np.array(
            [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], dtype=np.float64
        ),
    )


@pytest.mark.parametrize(
    ["r", "mu", "expected_distance", "expected_delta_shell"],
    [
        # Radially outward: the path is a radius, d = r_outer - r.
        (7.5e14, 1.0, SHELL_R_OUTER - 7.5e14, 1),
        # Radially inward: d = r - r_inner, and the packet enters the shell
        # below.
        (7.5e14, -1.0, 7.5e14 - SHELL_R_INNER, -1),
        # Perpendicular to the radius: right triangle with hypotenuse r_outer,
        # d = sqrt(r_outer^2 - r^2). Impact parameter r > r_inner, so outward.
        (7.5e14, 0.0, np.sqrt(SHELL_R_OUTER**2 - 7.5e14**2), 1),
        # Sitting on the outer boundary heading outward: zero distance.
        (SHELL_R_OUTER, 0.5, 0.0, 1),
        # Sitting on the inner boundary heading inward: zero distance.
        (SHELL_R_INNER, -0.5, 0.0, -1),
        # Starting on the outer boundary heading inward at mu = -0.5. The
        # impact parameter r_outer * sqrt(1 - 0.25) = 7.48e14 cm exceeds
        # r_inner, so the ray misses the inner sphere and re-exits through
        # r_outer along a chord of length 2 * r_outer * |mu| = r_outer.
        (SHELL_R_OUTER, -0.5, SHELL_R_OUTER, 1),
    ],
)
def test_calculate_distance_boundary_exact_geometry(
    r: float, mu: float, expected_distance: float, expected_delta_shell: int
) -> None:
    """
    Claim: the distance to the shell boundary matches closed-form geometry.

    Regime: straight-line propagation inside a single spherical shell.
    Verification: radial, perpendicular, on-boundary and chord cases whose
    distances follow from elementary geometry, independent of the quadratic
    roots used by the implementation.
    """
    distance, delta_shell = calculate_distances.calculate_distance_boundary(
        r, mu, SHELL_R_INNER, SHELL_R_OUTER
    )

    assert delta_shell == expected_delta_shell, (
        "Packet would be sent to the wrong neighbouring shell"
    )
    assert_allclose(
        distance, expected_distance, rtol=1e-14, atol=DISTANCE_ATOL_CM
    )


def test_calculate_distance_boundary_lands_on_target_sphere() -> None:
    """
    Claim: for any packet inside the shell, the returned distance is
    non-negative, the packet lands on the boundary sphere indicated by
    ``delta_shell`` after travelling it, and ``delta_shell`` is -1 exactly
    when the ray intersects the inner sphere.

    Regime: r_inner < r < r_outer, mu uniform in [-1, 1].
    Verification: the landing radius is recomputed with the law of cosines
    (the forward map the solver inverts), and the inner-sphere hit criterion
    is the impact-parameter test b = r * sqrt(1 - mu^2) < r_inner for an
    inward ray, rather than the discriminant form used by the implementation.
    An inward ray meets the inner sphere before its point of closest approach
    to the centre (at path length -r * mu), which rules out the far-side root.
    """
    # Arbitrary fixed seed so the 1000 sampled rays are deterministic.
    rng = np.random.default_rng(seed=2918)
    radii = rng.uniform(SHELL_R_INNER, SHELL_R_OUTER, size=1000)
    mus = rng.uniform(-1.0, 1.0, size=1000)

    for r, mu in zip(radii, mus, strict=True):
        distance, delta_shell = calculate_distances.calculate_distance_boundary(
            r, mu, SHELL_R_INNER, SHELL_R_OUTER
        )

        impact_parameter = r * np.sqrt(1.0 - mu * mu)
        hits_inner = mu < 0.0 and impact_parameter < SHELL_R_INNER
        expected_delta_shell = -1 if hits_inner else 1
        target_radius = SHELL_R_INNER if hits_inner else SHELL_R_OUTER
        r_landing, _ = propagate_along_ray(r, mu, distance)

        assert distance >= 0.0, f"Negative path length at r={r}, mu={mu}"
        assert delta_shell == expected_delta_shell, (
            f"Wrong shell crossing at r={r}, mu={mu}"
        )
        assert_allclose(
            r_landing,
            target_radius,
            rtol=LANDING_RTOL,
            err_msg=f"Packet does not land on the boundary at r={r}, mu={mu}",
        )
        if hits_inner:
            assert distance <= -r * mu, (
                f"Far-side intersection with the inner sphere at r={r}, mu={mu}"
            )


@pytest.mark.parametrize(
    ["mu_offset", "expected_delta_shell"],
    [
        # 1e-6 steeper than the tangent direction: the ray clips the inner
        # sphere.
        (-1e-6, -1),
        # 1e-6 shallower than the tangent direction: the ray passes the inner
        # sphere and exits through r_outer.
        (1e-6, 1),
    ],
)
def test_calculate_distance_boundary_grazing_inner_sphere(
    mu_offset: float, expected_delta_shell: int
) -> None:
    """
    Claim: rays just either side of the tangent to the inner sphere are sent
    to the correct neighbouring shell and land on that shell's boundary.

    Regime: near-grazing incidence, where the discriminant
    r_inner^2 - r^2 (1 - mu^2) passes through zero and is most sensitive to
    rounding.
    Verification: the tangent direction mu_t = -sqrt(1 - (r_inner / r)^2)
    gives the switch-over analytically; the landing radius is checked with
    the law of cosines.
    """
    r = 7.5e14
    mu_tangent = -np.sqrt(1.0 - (SHELL_R_INNER / r) ** 2)
    mu = mu_tangent + mu_offset

    distance, delta_shell = calculate_distances.calculate_distance_boundary(
        r, mu, SHELL_R_INNER, SHELL_R_OUTER
    )
    target_radius = (
        SHELL_R_INNER if expected_delta_shell == -1 else SHELL_R_OUTER
    )
    r_landing, _ = propagate_along_ray(r, mu, distance)

    assert delta_shell == expected_delta_shell
    assert_allclose(r_landing, target_radius, rtol=LANDING_RTOL)


@pytest.mark.xfail(
    strict=True,
    reason=(
        "mu == 0 takes the inward branch of calculate_distance_boundary; on "
        "r == r_inner the discriminant rounds slightly positive and returns a "
        "negative distance with delta_shell = -1."
    ),
)
@pytest.mark.parametrize("mu", [0.0, -0.0])
def test_calculate_distance_boundary_tangent_on_inner_boundary(
    mu: float,
) -> None:
    """
    Claim: a packet on the inner boundary moving tangentially (mu = 0) moves
    away from the inner sphere, so it must exit through r_outer after
    d = sqrt(r_outer^2 - r_inner^2).

    Regime: measure-zero edge case at the inner boundary; documents a known
    defect rather than a frequently reached state.
    Verification: right triangle with legs r_inner and d and hypotenuse
    r_outer.
    """
    distance, delta_shell = calculate_distances.calculate_distance_boundary(
        SHELL_R_INNER, mu, SHELL_R_INNER, SHELL_R_OUTER
    )

    assert delta_shell == 1
    assert_allclose(
        distance,
        np.sqrt(SHELL_R_OUTER**2 - SHELL_R_INNER**2),
        rtol=1e-14,
        atol=DISTANCE_ATOL_CM,
    )


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
    packet_params, expected_params, static_packet, time_explosion
):
    nu_line = packet_params["nu_line"]
    is_last_line = packet_params["is_last_line"]

    velocity = static_packet.r / time_explosion
    doppler_factor = frame_transformations.get_doppler_factor(
        velocity, static_packet.mu, False
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


@pytest.mark.parametrize("enable_full_relativity", [False, True])
# beta = v / c of the packet's starting position: 0.03 (9000 km/s) is typical
# of the photospheric region of a Type Ia supernova; 0.2 stresses the
# relativistic terms.
@pytest.mark.parametrize("beta", [0.03, 0.2])
@pytest.mark.parametrize("mu", [-1.0, -0.5, 0.0, 0.5, 1.0])
def test_calculate_distance_line_reaches_resonance(
    mu: float, beta: float, enable_full_relativity: bool
) -> None:
    """
    Claim: after travelling the returned line distance, the packet's
    comoving-frame frequency equals the line frequency.

    Regime: homologous expansion v = r / t, partial (first-order) and full
    special relativity, all propagation directions.
    Verification: the packet is advanced with the law of cosines and the
    comoving frequency at the new position is computed from the Doppler
    factor written out explicitly below, rather than from the distance
    formula under test.
    """
    time_explosion = 13.0 * 86400.0  # 13 days [s]
    nu_lab = 6.0e14  # ~500 nm, optical [Hz]
    ct = C_SPEED_OF_LIGHT * time_explosion
    r = beta * ct  # homologous: v = r / t

    def comoving_nu(r_now: float, mu_now: float) -> float:
        # nu_cmf = nu_lab * (1 - beta * mu) to first order, times the Lorentz
        # factor gamma = 1 / sqrt(1 - beta^2) in full relativity.
        beta_now = r_now / ct
        doppler_factor = 1.0 - beta_now * mu_now
        if enable_full_relativity:
            doppler_factor /= np.sqrt(1.0 - beta_now * beta_now)
        return nu_lab * doppler_factor

    comov_nu_start = comoving_nu(r, mu)
    # Place the line 0.1% redward of the current comoving frequency
    # (~300 km/s), far above CLOSE_LINE_THRESHOLD (1e-14), so the solver
    # takes the propagation branch rather than the zero-distance shortcut.
    nu_line = comov_nu_start * (1.0 - 1e-3)
    packet = radiative_packet.RPacket(
        r=r, mu=mu, nu=nu_lab, energy=1.0, seed=1963
    )

    distance = calculate_distances.calculate_distance_line(
        packet,
        comov_nu_start,
        False,
        nu_line,
        time_explosion,
        enable_full_relativity,
    )
    r_new, mu_new = propagate_along_ray(r, mu, distance)

    assert distance > 0.0
    # A rounding error delta in the distance shifts the comoving frequency by
    # ~nu * delta / ct. Even the cancellation in the full-relativity root
    # (losing ~3 digits for a 1e-3 offset) leaves this near 1e-16, so 1e-12
    # only fails on a wrong formula.
    assert_allclose(
        comoving_nu(r_new, mu_new),
        nu_line,
        rtol=1e-12,
        err_msg="Packet is not in resonance with the line after moving",
    )


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
    update_estimators_line(
        estimators,
        static_packet,
        cur_line_id,
        distance_trace,
        time_explosion,
        enable_full_relativity=False,
    )

    assert_allclose(estimators.mean_intensity_blueward, expected_j_blue)
    assert_allclose(estimators.energy_deposition_line_rate, expected_Edotlu)


# TODO set RNG consistently
# TODO: update this test to use the correct trace_packet
@pytest.mark.xfail(reason="Need to fix estimator differences across runs")
def test_trace_packet(
    packet,
    verysimple_time_explosion,
    verysimple_opacity_state,
    verysimple_geometry,
    verysimple_estimators_line,
    set_seed_fixture,
):
    set_seed_fixture(1963)
    packet.initialize_line_id(
        verysimple_opacity_state, verysimple_time_explosion
    )
    distance, interaction_type, delta_shell = r_packet_transport.trace_packet(
        packet,
        verysimple_geometry,
        verysimple_time_explosion,
        verysimple_opacity_state,
        verysimple_estimators_line,
        chi_continuum=1.0,  # Placeholder value
        escat_prob=0.5,  # Placeholder value
        enable_full_relativity=False,
        disable_line_scattering=False,
    )

    assert delta_shell == 1
    assert interaction_type == InteractionType.LINE
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
    geometry,
    time_explosion,
    estimators,
    ENABLE_FULL_RELATIVITY,
):
    distance = 1.0e13
    packet.nu = packet_params["nu"]
    packet.mu = packet_params["mu"]
    packet.energy = packet_params["energy"]
    packet.r = packet_params["r"]

    montecarlo_globals.ENABLE_FULL_RELATIVITY = ENABLE_FULL_RELATIVITY
    move_r_packet.recompile()  # This must be done as move_r_packet was jitted with ENABLE_FULL_RELATIVITY
    velocity = packet.r / time_explosion
    doppler_factor = frame_transformations.get_doppler_factor(
        velocity, packet.mu, ENABLE_FULL_RELATIVITY
    )

    move_r_packet(
        packet,
        distance,
        geometry,
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

    montecarlo_globals.ENABLE_FULL_RELATIVITY = False
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
    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)
    assert packet.status == radiative_packet.PacketStatus.EMITTED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -133, 132), (132, -133, 133), (132, -1e9, 133)],
)
def test_move_packet_across_shell_boundary_reabsorbed(
    packet, current_shell_id, delta_shell, no_of_shells
):
    packet.current_shell_id = current_shell_id
    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)
    assert packet.status == radiative_packet.PacketStatus.REABSORBED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -1, 199), (132, 0, 132), (132, 20, 154)],
)
def test_move_packet_across_shell_boundary_increment(
    packet, current_shell_id, delta_shell, no_of_shells
):
    packet.current_shell_id = current_shell_id
    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)
    assert packet.current_shell_id == current_shell_id + delta_shell
