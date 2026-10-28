from collections.abc import Callable
from decimal import Decimal, localcontext

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.macro_atom.macroatom_solver import LegacyMacroAtomSolver
from tardis.opacities.opacity_state_numba import OpacityStateNumba
from tardis.opacities.tau_sobolev import (
    calculate_sobolev_line_opacity as classic_calculate_sobolev_line_opacity,
)
from tardis.transport.geometry.calculate_distances import (
    calculate_comoving_frequency_nonhomologous,
    calculate_distance_boundary,
    calculate_distance_line_nonhomologous,
    calculate_projected_gradient_zero_distances,
    get_line_id_range_nonhomologous,
)
from tardis.transport.montecarlo.configuration.constants import (
    C_SPEED_OF_LIGHT,
    CLOSE_LINE_THRESHOLD,
)
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
)
from tardis.transport.montecarlo.modes.nonhomologous.interaction_event_callers import (
    line_scatter_event as nonhomologous_line_scatter_event,
)
from tardis.transport.montecarlo.modes.nonhomologous.interaction_events import (
    LineInteractionType,
)
from tardis.transport.montecarlo.modes.nonhomologous.interaction_events import (
    line_emission as nonhomologous_line_emission,
)
from tardis.transport.montecarlo.modes.nonhomologous.interaction_events import (
    thomson_scatter as nonhomologous_thomson_scatter,
)
from tardis.transport.montecarlo.modes.nonhomologous.opacity_solver import (
    OpacitySolver as NonhomologousOpacitySolver,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    trace_packet as nonhomologous_trace_packet,
)
from tardis.transport.montecarlo.modes.nonhomologous.tau_sobolev import (
    calculate_beta_sobolev,
    calculate_beta_sobolev_directional,
)
from tardis.transport.montecarlo.modes.nonhomologous.tau_sobolev import (
    calculate_sobolev_line_opacity as nonhomologous_calculate_sobolev_line_opacity,
)
from tardis.transport.montecarlo.nonhomologous_grid import (
    FLOAT64_ROUNDOFF_TOLERANCE,
    depressed_quartic,
    solve_resonance_quartic,
)
from tardis.transport.montecarlo.packets.movement import (
    move_packet_across_shell_boundary,
    move_r_packet,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
    RPacket,
)
from tardis.transport.montecarlo.utils import MonteCarloException


@pytest.mark.parametrize(
    ["A", "B", "C", "D", "E", "expected_roots"],
    [
        # x^4 - 14x^3 + 71x^2 - 154x + 120 = 0
        # roots 2, 3, 4, 5
        (1.0, -14.0, 71.0, -154.0, 120.0, [2.0, 3.0, 4.0, 5.0]),
        # x^4 - x^3 = 0
        # x^3(x - 1) = 0; only real root other than 0 is 1.0
        (1.0, -1.0, 0.0, 0.0, 0.0, [0.0, 0.0, 0.0, 1.0]),
        # x^4 - 10x^3 + 35x^2 - 50x + 24 = 0
        # roots 1, 2, 3, 4
        (1.0, -10.0, 35.0, -50.0, 24.0, [1.0, 2.0, 3.0, 4.0]),
        # x^4 - 1 = 0; two real and two complex roots
        (1.0, 0.0, 0.0, 0.0, -1.0, [-1.0, 1.0]),
        # x^4 - 5x^2 + 4 = 0; biquadratic with four real roots
        (1.0, 0.0, -5.0, 0.0, 4.0, [-2.0, -1.0, 1.0, 2.0]),
        # (x + 2)(x - 3)(x^2 + 1) = 0; two real roots
        (1.0, -1.0, -5.0, -1.0, -6.0, [-2.0, 3.0]),
        # (x^2 + 1)^2 = 0; repeated complex roots
        (1.0, 0.0, 2.0, 0.0, 1.0, []),
        # (x - 1)^2(x + 2)^2 = 0; two repeated real roots
        (1.0, 2.0, -3.0, -4.0, 4.0, [-2.0, -2.0, 1.0, 1.0]),
        # (x + 4)x^2(x - 1) = 0; repeated zero roots
        (1.0, 3.0, -4.0, 0.0, 0.0, [-4.0, 0.0, 0.0, 1.0]),
        # (x^2 - 1)(x^2 + 1e-14) = 0; do not make complex roots real
        (1.0, 0.0, -1.0 + 1.0e-14, 0.0, -1.0e-14, [-1.0, 1.0]),
    ],
)
def test_depressed_quartic(
    A: float,
    B: float,
    C: float,
    D: float,
    E: float,
    expected_roots: list[float],
) -> None:
    """
    Standalone unit test to check accurate calculation of expected roots.
    Not a regression test.
    """
    roots = np.asarray(depressed_quartic(A, B, C, D, E))
    real_roots = np.sort(roots[np.isfinite(roots)])

    npt.assert_allclose(real_roots, expected_roots, rtol=1.0e-11, atol=1.0e-12)
    assert np.count_nonzero(np.isnan(roots)) == 4 - len(expected_roots)
    for root in real_roots:
        numerator = abs(((((A * root) + B) * root + C) * root + D) * root + E)
        denominator = (
            abs(A) * abs(root) ** 4
            + abs(B) * abs(root) ** 3
            + abs(C) * abs(root) ** 2
            + abs(D) * abs(root)
            + abs(E)
        )
        scaled_residual = numerator / denominator if denominator else numerator
        assert scaled_residual < 1.0e-11


def test_depressed_quartic_preserves_close_resonance_roots() -> None:
    """Do not merge distinct roots close to zero projected velocity."""
    coefficients = (
        1.0,
        -0.1133408,
        6400.0032115342365,
        -1133.408,
        32.1153423616,
    )
    roots = np.asarray(depressed_quartic(*coefficients))
    real_roots = np.sort(roots[np.isfinite(roots)])

    npt.assert_allclose(
        real_roots,
        [0.035419000833124964, 0.14167578672128847],
        rtol=1.0e-11,
        atol=1.0e-13,
    )


def test_depressed_quartic_homologous_repeated_root() -> None:
    """Recover the repeated resonance root without a broad root clamp."""
    line_velocity = 1.0e-3
    impact_parameter = 200.0
    coefficients = (
        1.0,
        -2.0 * line_velocity,
        line_velocity**2 + impact_parameter**2,
        -2.0 * line_velocity * impact_parameter**2,
        line_velocity**2 * impact_parameter**2,
    )
    roots = np.asarray(depressed_quartic(*coefficients))
    real_roots = np.sort(roots[np.isfinite(roots)])

    npt.assert_allclose(
        real_roots,
        [line_velocity, line_velocity],
        rtol=1.0e-12,
        atol=1.0e-15,
    )


def test_depressed_quartic_does_not_force_complex_roots_to_zero() -> None:
    """Return no real roots for a positive biquadratic."""
    roots = np.asarray(depressed_quartic(1.0, 0.0, 5.0, 0.0, 4.0))

    assert np.all(np.isnan(roots))


@pytest.mark.parametrize(
    (
        "line_velocity",
        "impact_parameter_squared",
        "velocity_intercept",
        "expected_roots",
    ),
    [
        (1.0e-3, 40000.0, 0.0, [1.0e-3, 1.0e-3]),
        (0.0, 10000.0, -60.0, [0.0, 0.0]),
        (
            0.0,
            10000.0,
            -120.0,
            [-np.sqrt(4400.0), 0.0, 0.0, np.sqrt(4400.0)],
        ),
        (0.0, 10000.0, -100.0, [0.0, 0.0, 0.0, 0.0]),
        (0.25, 0.0, -0.75, [-0.5, 0.0, 0.0, 1.0]),
    ],
)
def test_solve_resonance_quartic_degenerate_cases(
    line_velocity: float,
    impact_parameter_squared: float,
    velocity_intercept: float,
    expected_roots: list[float],
) -> None:
    """Solve exact homologous, line-center, and radial limits directly."""
    roots = np.asarray(
        solve_resonance_quartic(
            line_velocity,
            impact_parameter_squared,
            velocity_intercept,
        )
    )
    real_roots = np.sort(roots[np.isfinite(roots)])

    npt.assert_allclose(real_roots, expected_roots)
    assert np.count_nonzero(np.isnan(roots)) == 4 - len(expected_roots)


def test_solve_resonance_quartic_preserves_physical_and_squared_roots() -> None:
    """Keep the physical and sign-reversed roots distinct near line center."""
    line_velocity = 0.0566704
    impact_parameter_squared = 10000.0
    velocity_intercept = -60.0
    roots = np.asarray(
        solve_resonance_quartic(
            line_velocity,
            impact_parameter_squared,
            velocity_intercept,
        )
    )
    real_roots = np.sort(roots[np.isfinite(roots)])

    npt.assert_allclose(
        real_roots,
        [0.035419000833124964, 0.14167578672128847],
        rtol=1.0e-11,
        atol=1.0e-13,
    )
    physical_residuals = (
        real_roots
        + velocity_intercept
        * real_roots
        / np.sqrt(real_roots**2 + impact_parameter_squared)
        - line_velocity
    )
    assert np.count_nonzero(np.abs(physical_residuals) < 1.0e-12) == 1


def test_solve_resonance_quartic_preserves_very_close_roots() -> None:
    """Resolve the stronger near-line-center case from the V838 model."""
    roots = np.asarray(solve_resonance_quartic(1.0e-3, 40000.0, -120.0))
    real_roots = np.sort(roots[np.isfinite(roots)])

    npt.assert_allclose(
        real_roots,
        [0.0006250000000011444, 0.0024999999997070313],
        rtol=1.0e-11,
        atol=1.0e-14,
    )


def test_solve_resonance_quartic_recovers_small_quadratic_factor() -> None:
    """Preserve near-zero roots when the velocity intercept dominates."""
    roots = np.asarray(
        solve_resonance_quartic(
            -2.529310910660824e-6,
            15.508081320848817**2,
            -210.42344701939626,
        )
    )

    npt.assert_allclose(
        np.sort(roots),
        [
            -209.85120332422667,
            -1.736134818207714e-7,
            2.012399543138253e-7,
            209.85119823797838,
        ],
        rtol=1.0e-11,
        atol=1.0e-14,
    )


def test_solve_resonance_quartic_recovers_tiny_shifted_roots() -> None:
    """Use reciprocal roots when the depressed-coordinate shift cancels."""
    roots = np.asarray(
        solve_resonance_quartic(
            -0.0113652295,
            0.00233329371**2,
            -973.326786,
        )
    )

    npt.assert_allclose(
        np.sort(roots),
        [
            -973.3381512267034,
            -2.72450685e-8,
            2.72451991e-8,
            973.3154207677032,
        ],
        rtol=1.0e-9,
        atol=1.0e-16,
    )


@pytest.mark.parametrize("coefficient_scale", [1.0e-100, 1.0, 1.0e100])
def test_depressed_quartic_coefficient_scale_invariance(
    coefficient_scale: float,
) -> None:
    """Normalize coefficients before evaluating powers and invariants."""
    coefficients = coefficient_scale * np.asarray(
        [1.0, -14.0, 71.0, -154.0, 120.0]
    )
    roots = np.asarray(depressed_quartic(*coefficients))

    npt.assert_allclose(
        np.sort(roots[np.isfinite(roots)]),
        [2.0, 3.0, 4.0, 5.0],
        rtol=1.0e-11,
        atol=1.0e-12,
    )


@pytest.mark.parametrize("root_scale", [1.0e-100, 1.0, 1.0e100])
def test_resonance_quartic_variable_scale_invariance(
    root_scale: float,
) -> None:
    """Scale resonance variables before evaluating their invariants."""
    roots = np.asarray(
        solve_resonance_quartic(
            0.0566704 * root_scale,
            10000.0 * root_scale**2,
            -60.0 * root_scale,
        )
    )
    real_roots = np.sort(roots[np.isfinite(roots)])

    npt.assert_allclose(
        real_roots,
        root_scale * np.asarray([0.035419000833124964, 0.14167578672128847]),
        rtol=1.0e-11,
        atol=1.0e-13 * root_scale,
    )


def test_nonhomologous_distance_solver_disables_fastmath() -> None:
    """Preserve finite checks and cancellation-sensitive root arithmetic."""
    assert depressed_quartic.targetoptions["fastmath"] is False
    assert solve_resonance_quartic.targetoptions["fastmath"] is False
    assert (
        calculate_distance_line_nonhomologous.targetoptions["fastmath"]
        is False
    )


def test_nonhomologous_distance_preserves_near_line_center_root() -> None:
    """Select the physical member of a close quartic root pair."""
    shell_width = 1.0e12
    geometry = NumbaRadial1DGeometry(
        np.asarray([1.0e14]),
        np.asarray([1.01e14]),
        np.asarray([1.0e7]),
        np.asarray([1.025e7]),
    )
    packet = RPacket(
        r=1.0e14,
        mu=0.0,
        nu=1.0e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0
    line_velocity = 0.0566704
    line_frequency = packet.nu * (
        1.0 - line_velocity * 2.5e5 / C_SPEED_OF_LIGHT
    )

    distance = calculate_distance_line_nonhomologous(
        packet, geometry, line_frequency
    )

    npt.assert_allclose(
        distance,
        0.14167578672128847 * shell_width,
        rtol=1.0e-9,
    )


def test_nonhomologous_distance_rejects_sign_reversed_root() -> None:
    """Select a root constructed from the unsquared resonance equation."""
    shell_width = 1.0e12
    shell_velocity_difference = 2.5e5
    geometry = NumbaRadial1DGeometry(
        np.asarray([2.0e14]),
        np.asarray([2.01e14]),
        np.asarray([2.0e7]),
        np.asarray([2.025e7]),
    )
    packet = RPacket(
        r=2.0e14,
        mu=0.0,
        nu=1.0e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0

    # Choose a physical projected position x, then construct the target line
    # velocity N from the original equation
    # N = x + Q*x/sqrt(x**2 + H**2). Squaring this equation to form the
    # quartic also produces a closer, sign-reversed root that must be rejected.
    expected_scaled_projected_position = 2.5e-3
    scaled_impact_parameter = packet.r / shell_width
    velocity_gradient = shell_velocity_difference / shell_width
    velocity_intercept = (
        geometry.v_outer[0] - velocity_gradient * geometry.r_outer[0]
    )
    scaled_velocity_intercept = (
        velocity_intercept / shell_velocity_difference
    )
    line_velocity = expected_scaled_projected_position + (
        scaled_velocity_intercept
        * expected_scaled_projected_position
        / np.sqrt(
            expected_scaled_projected_position**2
            + scaled_impact_parameter**2
        )
    )
    line_frequency = packet.nu * (
        1.0
        - line_velocity * shell_velocity_difference / C_SPEED_OF_LIGHT
    )

    distance = calculate_distance_line_nonhomologous(
        packet, geometry, line_frequency
    )

    npt.assert_allclose(
        distance,
        expected_scaled_projected_position * shell_width,
        rtol=1.0e-8,
    )


def test_nonhomologous_line_reappears_after_gradient_reversal(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
) -> None:
    """Verify two physical resonances on opposite sides of a gradient zero.

    Claim
    -----
    One line can be crossed twice when the projected velocity reverses.

    Regime
    ------
    Piecewise-linear radial velocity, nonradial ray, and negligible opacity.

    Verification
    ------------
    The unsquared projected-velocity equation fixes both resonance positions;
    its derivative independently fixes the two estimator contributions.
    """
    impact_parameter = 1.0e14  # cm
    first_x = 0.25 * impact_parameter
    second_x = 1.25 * impact_parameter
    velocity_intercept = 2.0e7  # cm / s
    first_radius = np.hypot(first_x, impact_parameter)
    second_radius = np.hypot(second_x, impact_parameter)
    velocity_gradient = velocity_intercept * (
        second_x / second_radius - first_x / first_radius
    ) / (first_x - second_x)

    inner_radius = 0.9 * impact_parameter
    outer_radius = 1.7 * impact_parameter
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([velocity_gradient * inner_radius + velocity_intercept]),
        np.asarray([velocity_gradient * outer_radius + velocity_intercept]),
    )
    packet = RPacket(
        r=impact_parameter,
        mu=0.0,
        nu=1.0e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0
    first_projected_velocity = (
        velocity_gradient * first_x
        + velocity_intercept * first_x / first_radius
    )
    line_frequency = packet.nu * (
        1.0 - first_projected_velocity / C_SPEED_OF_LIGHT
    )
    classic_opacity_state.line_list_nu = np.asarray([line_frequency, 1.0e14])

    boundary_distance, _ = calculate_distance_boundary(
        packet.r, packet.mu, inner_radius, outer_radius
    )
    first_zero, _, zero_count = calculate_projected_gradient_zero_distances(
        packet, geometry, boundary_distance
    )
    turning_radius = (
        -velocity_intercept * impact_parameter**2 / velocity_gradient
    ) ** (1.0 / 3.0)
    expected_zero = np.sqrt(turning_radius**2 - impact_parameter**2)
    assert zero_count == 1
    npt.assert_allclose(first_zero, expected_zero, rtol=0.0, atol=100.0)

    first_distance = calculate_distance_line_nonhomologous(
        packet, geometry, line_frequency, 0.0, first_zero
    )
    second_distance = calculate_distance_line_nonhomologous(
        packet, geometry, line_frequency, first_zero, boundary_distance
    )
    # A frequency ULP near 1e15 Hz maps to less than 100 cm at these gradients.
    npt.assert_allclose(
        [first_distance, second_distance],
        [first_x, second_x],
        rtol=0.0,
        atol=100.0,
    )

    _, interaction_type, _, _, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        True,
        0,
        0,
    )
    assert interaction_type == InteractionType.BOUNDARY

    first_gradient = velocity_gradient + (
        velocity_intercept * impact_parameter**2 / first_radius**3
    )
    second_gradient = velocity_gradient + (
        velocity_intercept * impact_parameter**2 / second_radius**3
    )
    assert first_gradient > 0.0 > second_gradient
    expected_energy = line_frequency / packet.nu
    expected_deposition = 2.0 * expected_energy
    expected_blueward = (
        expected_energy
        / packet.nu
        * (1.0 / abs(first_gradient) + 1.0 / abs(second_gradient))
    )
    npt.assert_allclose(
        line_estimators.energy_deposition_line_rate[0, 0],
        expected_deposition,
        rtol=1.0e-10,
    )
    npt.assert_allclose(
        line_estimators.mean_intensity_blueward[0, 0],
        expected_blueward,
        rtol=1.0e-10,
    )


@pytest.mark.parametrize("frequency_ulp_shift", [0, 1])
def test_nonhomologous_line_at_projected_gradient_zero_is_singular(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
    frequency_ulp_shift: int,
) -> None:
    """A line within rounding uncertainty of a turning point is unresolved.

    The quartic can split a tangent resonance into two apparent crossings, one
    on either side of the analytic projected-gradient zero. The line search
    must identify the represented extremum before accepting either root. A
    frequency one ULP away is also unresolved at this ill-conditioned point.
    """
    impact_parameter = 1.0e14  # cm
    first_x = 0.25 * impact_parameter
    second_x = 1.25 * impact_parameter
    velocity_intercept = 2.0e7  # cm / s
    first_radius = np.hypot(first_x, impact_parameter)
    second_radius = np.hypot(second_x, impact_parameter)
    velocity_gradient = velocity_intercept * (
        second_x / second_radius - first_x / first_radius
    ) / (first_x - second_x)
    inner_radius = 0.9 * impact_parameter
    outer_radius = 1.7 * impact_parameter
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([velocity_gradient * inner_radius + velocity_intercept]),
        np.asarray([velocity_gradient * outer_radius + velocity_intercept]),
    )
    packet = RPacket(
        r=impact_parameter,
        mu=0.0,
        nu=1.0e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0
    boundary_distance, _ = calculate_distance_boundary(
        packet.r, packet.mu, inner_radius, outer_radius
    )
    zero_distance, _, zero_count = calculate_projected_gradient_zero_distances(
        packet, geometry, boundary_distance
    )
    assert zero_count == 1
    turning_frequency = calculate_comoving_frequency_nonhomologous(
        packet, geometry, zero_distance
    )
    if frequency_ulp_shift:
        turning_frequency = np.nextafter(turning_frequency, np.inf)
    classic_opacity_state.line_list_nu = np.asarray(
        [turning_frequency, 1.0e14]
    )

    with pytest.raises(
        MonteCarloException,
        match="Sobolev optical depth is unresolved",
    ):
        nonhomologous_trace_packet(
            packet,
            geometry,
            classic_opacity_state,
            line_estimators,
            1.0e-30,
            False,
            True,
            0,
            0,
        )

    # A Thomson event before the extremum still takes precedence.
    distance, interaction_type, _, _, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-3,
        False,
        True,
        0,
        0,
    )
    assert interaction_type == InteractionType.ESCATTERING
    assert distance < zero_distance


def test_nonhomologous_line_below_projected_frequency_minimum_is_bounded(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
) -> None:
    """A below-minimum line is unresolved nearby and missed farther away.

    The represented velocity field has a positive-to-negative projected
    gradient reversal, so its comoving frequency has a minimum. Decimal
    arithmetic independently places the chosen line below that minimum. The
    near line is within the float64 conditioning band; the farther line is
    outside it and cannot contribute an estimator.
    """
    impact_parameter = 1.0e14  # cm
    first_x = 0.25 * impact_parameter
    second_x = 1.25 * impact_parameter
    velocity_intercept = 2.0e7  # cm / s
    first_radius = np.hypot(first_x, impact_parameter)
    second_radius = np.hypot(second_x, impact_parameter)
    velocity_gradient = velocity_intercept * (
        second_x / second_radius - first_x / first_radius
    ) / (first_x - second_x)
    inner_radius = 0.9 * impact_parameter
    outer_radius = 1.7 * impact_parameter
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([velocity_gradient * inner_radius + velocity_intercept]),
        np.asarray([velocity_gradient * outer_radius + velocity_intercept]),
    )
    packet = RPacket(
        r=impact_parameter,
        mu=0.0,
        nu=1.0e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0
    boundary_distance, _ = calculate_distance_boundary(
        packet.r, packet.mu, inner_radius, outer_radius
    )
    zero_distance, _, zero_count = calculate_projected_gradient_zero_distances(
        packet, geometry, boundary_distance
    )
    assert zero_count == 1
    turning_frequency = calculate_comoving_frequency_nonhomologous(
        packet, geometry, zero_distance
    )
    unreachable_frequency = np.nextafter(turning_frequency, -np.inf)

    with localcontext() as decimal_context:
        decimal_context.prec = 60
        radius = Decimal(float(impact_parameter))
        gradient = Decimal(float(geometry.velocity_gradient[0]))
        intercept = Decimal(float(geometry.v_inner[0])) - gradient * Decimal(
            float(geometry.r_inner[0])
        )
        turning_radius = (
            -intercept * radius**2 / gradient
        ) ** (Decimal(1) / Decimal(3))
        turning_x = (turning_radius**2 - radius**2).sqrt()
        minimum_frequency = Decimal(float(packet.nu)) * (
            Decimal(1)
            - (
                gradient * turning_x
                + intercept * turning_x / turning_radius
            )
            / Decimal(float(C_SPEED_OF_LIGHT))
        )
        assert Decimal(float(unreachable_frequency)) < minimum_frequency

    classic_opacity_state.line_list_nu = np.asarray(
        [unreachable_frequency, 1.0e14]
    )
    with pytest.raises(
        MonteCarloException,
        match="Sobolev optical depth is unresolved",
    ):
        nonhomologous_trace_packet(
            packet,
            geometry,
            classic_opacity_state,
            line_estimators,
            1.0e-30,
            False,
            True,
            0,
            0,
        )
    assert line_estimators.energy_deposition_line_rate[0, 0] == 0.0
    assert line_estimators.mean_intensity_blueward[0, 0] == 0.0

    classic_opacity_state.line_list_nu[0] = turning_frequency - 8.0
    _, interaction_type, _, _, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        True,
        0,
        0,
    )
    assert interaction_type == InteractionType.BOUNDARY
    assert line_estimators.energy_deposition_line_rate[0, 0] == 0.0
    assert line_estimators.mean_intensity_blueward[0, 0] == 0.0


def test_nonhomologous_search_resolves_two_close_lines(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
) -> None:
    """Keep two distinct nearby lines in a monotonic projected flow.

    Claim
    -----
    The line search counts each of two resolvably close resonances once.

    Regime
    ------
    Positive projected gradient in a genuinely nonhomologous shell.

    Verification
    ------------
    Construct both rest frequencies from the unsquared ray equation at known
    distances; their local Doppler factors set independent energy estimates.
    """
    inner_radius = 1.0e14  # cm
    outer_radius = 2.0e14  # cm
    velocity_gradient = 8.0e-7  # s^-1
    velocity_intercept = 2.0e7  # cm / s
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([velocity_gradient * inner_radius + velocity_intercept]),
        np.asarray([velocity_gradient * outer_radius + velocity_intercept]),
    )
    packet = RPacket(
        r=1.25e14,
        mu=0.1,
        nu=1.0e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0
    first_distance = 1.0e12  # cm
    second_distance = first_distance + 1.0e7  # cm
    initial_x = packet.r * packet.mu
    impact_parameter_squared = packet.r**2 * (1.0 - packet.mu**2)
    x_positions = initial_x + np.asarray([first_distance, second_distance])
    resonance_radii = np.sqrt(x_positions**2 + impact_parameter_squared)
    projected_velocities = (
        velocity_gradient * x_positions
        + velocity_intercept * x_positions / resonance_radii
    )
    line_frequencies = packet.nu * (
        1.0 - projected_velocities / C_SPEED_OF_LIGHT
    )
    assert line_frequencies[0] > line_frequencies[1]
    classic_opacity_state.line_list_nu = line_frequencies

    solved_distances = np.asarray(
        [
            calculate_distance_line_nonhomologous(
                packet, geometry, line_frequency
            )
            for line_frequency in line_frequencies
        ]
    )
    npt.assert_allclose(
        solved_distances,
        [first_distance, second_distance],
        rtol=0.0,
        atol=100.0,
    )

    _, interaction_type, _, _, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        True,
        0,
        0,
    )
    assert interaction_type == InteractionType.BOUNDARY
    npt.assert_allclose(
        line_estimators.energy_deposition_line_rate[:, 0],
        line_frequencies / packet.nu,
        rtol=1.0e-10,
    )


def test_nonhomologous_thomson_scatter_keeps_same_direction_line_frontier(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
    bulk_estimators: EstimatorsBulk,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """Carry a crossed line through an elastic scattering on a monotonic ray.

    Claim
    -----
    Thomson scattering preserves the comoving frequency at the event and the
    next lower-frequency line remains reachable without revisiting the first.

    Regime
    ------
    Positive projected gradient for both sampled directions in one shell.

    Verification
    ------------
    Use analytic projected velocities to place one line before scattering and
    the second line a known distance after it; verify both estimator channels.
    """
    inner_radius = 1.0e14  # cm
    outer_radius = 2.0e14  # cm
    velocity_gradient = 8.0e-7  # s^-1
    velocity_intercept = 2.0e7  # cm / s
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([velocity_gradient * inner_radius + velocity_intercept]),
        np.asarray([velocity_gradient * outer_radius + velocity_intercept]),
    )
    packet = RPacket(r=1.25e14, mu=0.1, nu=1.0e15, energy=1.0, seed=42)
    packet.current_shell_id = 0
    first_line_distance = 1.0e12  # cm
    electron_distance = 2.0e12  # cm
    next_line_distance = 1.0e12  # cm after Thomson scattering

    initial_x = packet.r * packet.mu
    initial_b_squared = packet.r**2 * (1.0 - packet.mu**2)
    first_line_x = initial_x + first_line_distance
    first_line_radius = np.sqrt(first_line_x**2 + initial_b_squared)
    first_line_projected_velocity = (
        velocity_gradient * first_line_x
        + velocity_intercept * first_line_x / first_line_radius
    )
    first_line_frequency = packet.nu * (
        1.0 - first_line_projected_velocity / C_SPEED_OF_LIGHT
    )

    electron_x = initial_x + electron_distance
    electron_radius = np.sqrt(electron_x**2 + initial_b_squared)
    electron_mu = electron_x / electron_radius
    electron_velocity = velocity_gradient * electron_radius + velocity_intercept
    comov_frequency_at_electron = packet.nu * (
        1.0 - electron_velocity * electron_mu / C_SPEED_OF_LIGHT
    )
    seed = 42
    first_random, direction_random = np.random.RandomState(seed).random_sample(2)
    sampled_outgoing_mu = 2.0 * direction_random - 1.0
    outgoing_lab_frequency = comov_frequency_at_electron / (
        1.0 - electron_velocity * sampled_outgoing_mu / C_SPEED_OF_LIGHT
    )
    outgoing_x = electron_radius * sampled_outgoing_mu + next_line_distance
    outgoing_b_squared = electron_radius**2 * (1.0 - sampled_outgoing_mu**2)
    next_line_radius = np.sqrt(outgoing_x**2 + outgoing_b_squared)
    next_line_projected_velocity = (
        velocity_gradient * outgoing_x
        + velocity_intercept * outgoing_x / next_line_radius
    )
    next_line_frequency = outgoing_lab_frequency * (
        1.0 - next_line_projected_velocity / C_SPEED_OF_LIGHT
    )
    assert first_line_frequency > next_line_frequency
    classic_opacity_state.line_list_nu = np.asarray(
        [first_line_frequency, next_line_frequency]
    )
    classic_opacity_state.tau_sobolev[1, 0] = 1.0e6

    set_seed_fixture(seed)
    electron_opacity = -np.log(first_random) / electron_distance
    (
        traced_distance,
        first_event,
        _,
        next_untraversed_line_id,
        line_search_direction,
    ) = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        electron_opacity,
        False,
        False,
        0,
        0,
    )
    assert first_event == InteractionType.ESCATTERING
    npt.assert_allclose(traced_distance, electron_distance, rtol=1.0e-12)
    assert (next_untraversed_line_id, line_search_direction) == (1, 1)
    first_line_energy = line_estimators.energy_deposition_line_rate[0, 0]
    npt.assert_allclose(
        first_line_energy,
        first_line_frequency / packet.nu,
        rtol=1.0e-10,
    )

    move_r_packet(packet, traced_distance, geometry, bulk_estimators, False)
    incoming_mu = packet.mu
    incoming_lab_frequency = packet.nu
    nonhomologous_thomson_scatter(packet, geometry, False)
    npt.assert_allclose(packet.mu, sampled_outgoing_mu, rtol=0.0, atol=1.0e-15)
    npt.assert_allclose(
        packet.nu * (1.0 - electron_velocity * packet.mu / C_SPEED_OF_LIGHT),
        incoming_lab_frequency
        * (1.0 - electron_velocity * incoming_mu / C_SPEED_OF_LIGHT),
        rtol=1.0e-14,
    )

    (
        second_distance,
        second_event,
        _,
        next_untraversed_line_id,
        line_search_direction,
    ) = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        False,
        next_untraversed_line_id,
        line_search_direction,
    )
    assert second_event == InteractionType.LINE
    assert packet.next_line_id == 1
    assert (next_untraversed_line_id, line_search_direction) == (2, 1)
    npt.assert_allclose(second_distance, next_line_distance, rtol=0.0, atol=100.0)
    npt.assert_allclose(
        line_estimators.energy_deposition_line_rate[:, 0],
        [first_line_energy, packet.energy * next_line_frequency / packet.nu],
        rtol=1.0e-10,
    )


def test_nonhomologous_frequency_changing_emission_starts_new_line_history(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
    bulk_estimators: EstimatorsBulk,
) -> None:
    """Allow a crossed line after emission at a higher rest frequency.

    Claim
    -----
    A new emitted frequency begins a new line traversal history, so a line
    crossed before the emission can be crossed again later on the same ray.

    Regime
    ------
    Positive projected gradient and a frequency-changing line emission.

    Verification
    ------------
    The unsquared ray equation sets both resonances at prescribed distances;
    the packet is moved and re-emitted with the production routines.
    """
    inner_radius = 1.0e14  # cm
    outer_radius = 2.0e14  # cm
    velocity_gradient = 8.0e-7  # s^-1
    velocity_intercept = 2.0e7  # cm / s
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([velocity_gradient * inner_radius + velocity_intercept]),
        np.asarray([velocity_gradient * outer_radius + velocity_intercept]),
    )
    initial_radius = 1.25e14  # cm
    initial_mu = 0.1
    absorbed_frequency = 1.0e15  # Hz
    first_line_distance = 1.0e12  # cm
    repeated_line_distance = 3.0e12  # cm after emission
    initial_x = initial_radius * initial_mu
    impact_parameter_squared = initial_radius**2 * (1.0 - initial_mu**2)
    first_x = initial_x + first_line_distance
    first_radius = np.sqrt(first_x**2 + impact_parameter_squared)
    first_mu = first_x / first_radius
    first_velocity = velocity_gradient * first_radius + velocity_intercept
    incoming_lab_frequency = absorbed_frequency / (
        1.0 - first_velocity * first_mu / C_SPEED_OF_LIGHT
    )
    packet = RPacket(
        r=initial_radius,
        mu=initial_mu,
        nu=incoming_lab_frequency,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0

    repeated_x = first_x + repeated_line_distance
    repeated_radius = np.sqrt(repeated_x**2 + impact_parameter_squared)
    repeated_mu = repeated_x / repeated_radius
    repeated_velocity = velocity_gradient * repeated_radius + velocity_intercept
    emitted_frequency = absorbed_frequency * (
        1.0 - first_velocity * first_mu / C_SPEED_OF_LIGHT
    ) / (1.0 - repeated_velocity * repeated_mu / C_SPEED_OF_LIGHT)
    initial_velocity = velocity_gradient * initial_radius + velocity_intercept
    initial_comov_frequency = incoming_lab_frequency * (
        1.0 - initial_velocity * initial_mu / C_SPEED_OF_LIGHT
    )
    assert emitted_frequency > initial_comov_frequency > absorbed_frequency
    classic_opacity_state.line_list_nu = np.asarray(
        [emitted_frequency, absorbed_frequency]
    )
    classic_opacity_state.tau_sobolev[1, 0] = 1.0e6

    (
        first_distance,
        first_event,
        _,
        next_untraversed_line_id,
        line_search_direction,
    ) = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        False,
        0,
        0,
    )
    assert first_event == InteractionType.LINE
    assert packet.next_line_id == 1
    npt.assert_allclose(first_distance, first_line_distance, rtol=0.0, atol=100.0)
    assert (next_untraversed_line_id, line_search_direction) == (2, 1)
    first_deposition = line_estimators.energy_deposition_line_rate[1, 0]

    move_r_packet(packet, first_distance, geometry, bulk_estimators, False)
    nonhomologous_line_emission(
        packet,
        0,
        geometry,
        classic_opacity_state,
        False,
    )
    # A new line is emitted at the packet position; it is already traversed.
    next_untraversed_line_id = 1
    line_search_direction = 1

    (
        second_distance,
        second_event,
        _,
        next_untraversed_line_id,
        line_search_direction,
    ) = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        False,
        next_untraversed_line_id,
        line_search_direction,
    )
    assert second_event == InteractionType.LINE
    assert packet.next_line_id == 1
    assert (next_untraversed_line_id, line_search_direction) == (2, 1)
    npt.assert_allclose(
        second_distance, repeated_line_distance, rtol=0.0, atol=100.0
    )
    npt.assert_allclose(
        line_estimators.energy_deposition_line_rate[1, 0] - first_deposition,
        packet.energy * absorbed_frequency / packet.nu,
        rtol=1.0e-10,
    )
    assert line_estimators.energy_deposition_line_rate[0, 0] == 0.0


@pytest.mark.parametrize(
    ("projected_gradient", "expected_line_ids"),
    [
        (1.0e-6, (0, 1, 2)),
        (-1.0e-6, (2, 1, 0)),
    ],
)
@pytest.mark.parametrize(
    "end_frequency",
    [1.0e15, np.nextafter(1.0e15, np.inf)],
)
def test_nonhomologous_range_uses_projected_gradient_direction(
    projected_gradient: float,
    expected_line_ids: tuple[int, int, int],
    end_frequency: float,
) -> None:
    """Order candidate lines by the known gradient when endpoints round close.

    Claim
    -----
    The sign of the projected gradient fixes line-list traversal direction.

    Regime
    ------
    Equal or one-ULP-separated endpoint frequencies near 1e15 Hz.

    Verification
    ------------
    A descending line list has known integer order for each gradient sign.
    """
    start_frequency = 1.0e15  # Hz
    line_frequencies = np.asarray(
        [start_frequency + 5.0, start_frequency, start_frequency - 5.0]
    )
    line_id_step = 1 if projected_gradient > 0.0 else -1

    start_line_id, stop_line_id, returned_step = get_line_id_range_nonhomologous(
        line_frequencies,
        start_frequency,
        end_frequency,
        line_id_step,
    )

    assert returned_step == line_id_step
    assert tuple(range(start_line_id, stop_line_id, returned_step)) == (
        expected_line_ids
    )


def test_nonhomologous_zero_finder_keeps_representable_near_start_reversal() -> None:
    """Keep a sign-changing zero even inside the shell-width distance cutoff.

    Claim
    -----
    A representable projected-gradient reversal partitions the ray.

    Regime
    ------
    Positive radial velocity and a zero about 5 cm from the packet start,
    inside a 10 cm shell-width cutoff. This test does not assert that a line
    frequency can resolve the very small spectral excursion there.

    Verification
    ------------
    Reconstruct the linear field from its edge values and solve the analytic
    cubic with ``numpy.cbrt``; check gradient signs on either side.
    """
    inner_radius = 1.0e14  # cm
    outer_radius = 1.1e15  # cm
    packet_radius = 1.25e14  # cm
    packet_direction_cosine = 0.95
    velocity_intercept = 2.0e7  # cm / s
    projected_start = packet_radius * packet_direction_cosine
    impact_parameter_squared = packet_radius**2 * (
        1.0 - packet_direction_cosine**2
    )
    target_zero_distance = 5.0  # cm
    target_turning_x = projected_start + target_zero_distance
    target_turning_radius = np.sqrt(
        target_turning_x**2 + impact_parameter_squared
    )
    velocity_gradient = (
        -velocity_intercept
        * impact_parameter_squared
        / target_turning_radius**3
    )
    inner_velocity = velocity_gradient * inner_radius + velocity_intercept
    outer_velocity = velocity_gradient * outer_radius + velocity_intercept
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([inner_velocity]),
        np.asarray([outer_velocity]),
    )
    packet = RPacket(
        r=packet_radius,
        mu=packet_direction_cosine,
        nu=1.0e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0

    represented_gradient = (outer_velocity - inner_velocity) / (
        outer_radius - inner_radius
    )
    represented_intercept = inner_velocity - represented_gradient * inner_radius
    turning_radius = np.cbrt(
        -represented_intercept * impact_parameter_squared / represented_gradient
    )
    expected_zero_distance = (
        np.sqrt(turning_radius**2 - impact_parameter_squared)
        - projected_start
    )
    distance_cutoff = CLOSE_LINE_THRESHOLD * (outer_radius - inner_radius)
    assert 100.0 * np.spacing(projected_start) < expected_zero_distance
    assert expected_zero_distance < distance_cutoff

    radius_after_zero = np.sqrt(
        (projected_start + 20.0) ** 2 + impact_parameter_squared
    )
    gradient_at_start = represented_gradient + (
        represented_intercept * impact_parameter_squared / packet_radius**3
    )
    gradient_after_zero = represented_gradient + (
        represented_intercept * impact_parameter_squared / radius_after_zero**3
    )
    assert gradient_at_start > 0.0 > gradient_after_zero

    boundary_distance, _ = calculate_distance_boundary(
        packet.r, packet.mu, inner_radius, outer_radius
    )
    first_zero, _, zero_count = calculate_projected_gradient_zero_distances(
        packet, geometry, boundary_distance
    )
    assert zero_count == 1
    npt.assert_allclose(first_zero, expected_zero_distance, rtol=0.0, atol=1.0)


def test_nonhomologous_zero_finder_resolves_near_tangent_pair() -> None:
    """Resolve two ray-coordinate zeros despite a small radial separation.

    Claim
    -----
    A sign-changing pair near the impact parameter has two forward zeros.

    Regime
    ------
    The turning radius lies about 0.5 cm above a 1e14 cm impact parameter,
    below the old 1 cm radial cutoff, but the zeros are about 2e7 cm apart
    along the ray. Their spectral excursion at 1e15 Hz is below 1e-9 Hz,
    far below the 0.125 Hz frequency spacing; this is a geometric test only.

    Verification
    ------------
    Reconstruct the represented linear field from its edge values and solve
    its cubic with ``numpy.cbrt``. One radius ULP near 1e14 cm is 0.015625 cm;
    its amplification by ``r/x`` near the zeros motivates 2e5 cm tolerance.
    """
    impact_parameter = 1.0e14  # cm
    inner_radius = impact_parameter - 1.0e8  # cm
    outer_radius = 2.0 * impact_parameter  # cm
    velocity_intercept = -2.0e7  # cm / s
    target_turning_radius = impact_parameter + 0.5  # cm
    velocity_gradient = (
        -velocity_intercept
        * impact_parameter**2
        / target_turning_radius**3
    )
    inner_velocity = velocity_gradient * inner_radius + velocity_intercept
    outer_velocity = velocity_gradient * outer_radius + velocity_intercept
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([inner_velocity]),
        np.asarray([outer_velocity]),
    )
    initial_x = -2.0e7  # cm
    packet_radius = np.sqrt(impact_parameter**2 + initial_x**2)
    packet = RPacket(
        r=packet_radius,
        mu=initial_x / packet_radius,
        nu=1.0e15,
        energy=1.0,
        seed=1963,
    )
    packet.current_shell_id = 0

    represented_gradient = (outer_velocity - inner_velocity) / (
        outer_radius - inner_radius
    )
    represented_intercept = outer_velocity - represented_gradient * outer_radius
    represented_impact_squared = packet.r**2 * (1.0 - packet.mu**2)
    turning_radius = np.cbrt(
        -represented_intercept * represented_impact_squared
        / represented_gradient
    )
    turning_x = np.sqrt(turning_radius**2 - represented_impact_squared)
    projected_start = packet.r * packet.mu
    expected_zeros = np.asarray(
        [-turning_x - projected_start, turning_x - projected_start]
    )
    radial_cutoff = CLOSE_LINE_THRESHOLD * (outer_radius - inner_radius)
    assert turning_radius - np.sqrt(represented_impact_squared) < radial_cutoff
    assert expected_zeros[1] - expected_zeros[0] > 1.0e7

    sample_x = np.asarray([initial_x, 0.0, -initial_x])
    sample_radius = np.sqrt(sample_x**2 + represented_impact_squared)
    sample_gradient = represented_gradient + (
        represented_intercept * represented_impact_squared / sample_radius**3
    )
    assert sample_gradient[0] > 0.0
    assert sample_gradient[1] < 0.0
    assert sample_gradient[2] > 0.0

    boundary_distance, _ = calculate_distance_boundary(
        packet.r, packet.mu, inner_radius, outer_radius
    )
    first_zero, second_zero, zero_count = (
        calculate_projected_gradient_zero_distances(
            packet, geometry, boundary_distance
        )
    )
    assert zero_count == 2
    npt.assert_allclose(
        [first_zero, second_zero], expected_zeros, rtol=0.0, atol=2.0e5
    )


def test_nonhomologous_calculate_beta_sobolevs(
    nb_simulation_verysimple, regression_data
):
    """
    Analogous to
    tardis/opacities/tests/test_tau_sobolev.py@test_calculate_beta_sobolevs
    """
    legacy_plasma = nb_simulation_verysimple.plasma

    # Testing only the nonhomologous beta sobolev calculation, so start with
    # tau sobolevs from classic mode
    tau_sobolevs = classic_calculate_sobolev_line_opacity(
        legacy_plasma.lines,
        legacy_plasma.level_number_density,
        legacy_plasma.time_explosion,
        legacy_plasma.stimulated_emission_factor,
    )

    actual = calculate_beta_sobolev(tau_sobolevs)
    expected = regression_data.sync_ndarray(actual)
    npt.assert_allclose(actual, expected)


def test_directional_beta_sobolev_splits_projected_gradient_zero() -> None:
    """Resolve the escape-probability cusp at an interior gradient zero."""
    optical_depth_coefficient = pd.DataFrame([[1.0]])
    velocity_gradient = np.asarray([-1.0]) / u.s
    velocity_over_radius = np.asarray([1.0]) / u.s

    actual = calculate_beta_sobolev_directional(
        optical_depth_coefficient,
        velocity_gradient,
        velocity_over_radius,
        quadrature_order=20,
    )

    npt.assert_allclose(actual.to_numpy(), [[0.44819011846544604]], rtol=1.0e-8)


def test_directional_beta_sobolev_matches_homologous_limit() -> None:
    """Recover angle-independent escape probabilities under homology."""
    optical_depth_coefficient = pd.DataFrame([[0.0, 0.5], [1.0, 10.0]])
    homologous_gradient = np.asarray([2.0, 4.0]) / u.s
    expected = calculate_beta_sobolev(
        optical_depth_coefficient / np.asarray([2.0, 4.0])
    )

    actual = calculate_beta_sobolev_directional(
        optical_depth_coefficient,
        homologous_gradient,
        homologous_gradient,
        quadrature_order=2,
    )

    npt.assert_allclose(actual.to_numpy(), expected.to_numpy())


def test_directional_beta_sobolev_handles_zero_velocity_gradient() -> None:
    """Return zero escape for opaque lines in a static velocity field."""
    optical_depth_coefficient = pd.DataFrame([[0.0], [1.0]])
    zero_gradient = np.asarray([0.0]) / u.s

    actual = calculate_beta_sobolev_directional(
        optical_depth_coefficient,
        zero_gradient,
        zero_gradient,
    )

    npt.assert_allclose(
        actual.to_numpy(),
        [[1.0], [0.0]],
        rtol=FLOAT64_ROUNDOFF_TOLERANCE,
        atol=0.0,
    )


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, 11, 132), (132, 1, 133), (132, 2, 133)],
)
def test_nonhomologous_move_packet_across_shell_boundary_emitted(
    current_shell_id, delta_shell, no_of_shells
):
    """
    Analogous to
    tardis/transport/montecarlo/tests/test_montecarlo.py@test_move_packet_across_shell_boundary_emitted
    """
    packet = RPacket(r=7.5e14, mu=0.3, nu=0.4, energy=0.9, seed=1963)
    packet.current_shell_id = current_shell_id
    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)
    assert packet.status == PacketStatus.EMITTED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -133, 132), (132, -133, 133), (132, -1e9, 133)],
)
def test_nonhomologous_move_packet_across_shell_boundary_reabsorbed(
    current_shell_id, delta_shell, no_of_shells
):
    """
    Analogous to
    tardis/transport/montecarlo/tests/test_montecarlo.py@test_move_packet_across_shell_boundary_reabsorbed
    """
    packet = RPacket(r=7.5e14, mu=0.3, nu=0.4, energy=0.9, seed=1963)
    packet.current_shell_id = current_shell_id
    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)
    assert packet.status == PacketStatus.REABSORBED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -1, 199), (132, 0, 132), (132, 20, 154)],
)
def test_nonhomologous_move_packet_across_shell_boundary_increment(
    current_shell_id, delta_shell, no_of_shells
):
    """
    Analogous to
    tardis/transport/montecarlo/tests/test_montecarlo.py@test_move_packet_across_shell_boundary_increment
    """
    packet = RPacket(r=7.5e14, mu=0.3, nu=0.4, energy=0.9, seed=1963)
    packet.current_shell_id = current_shell_id
    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)
    assert packet.current_shell_id == current_shell_id + delta_shell


def test_nonhomologous_thomson_scatter(
    packet, verysimple_numba_active_radial_1d_geometry
):
    """
    Analogous to
    tardis/transport/montecarlo/tests/test_interaction.py@test_thomson_scatter
    """
    init_mu = packet.mu
    init_nu = packet.nu
    init_energy = packet.energy

    nonhomologous_thomson_scatter(
        packet, verysimple_numba_active_radial_1d_geometry, False
    )

    assert np.abs(packet.mu - init_mu) > 1e-7
    assert np.abs(packet.nu - init_nu) > 1e-7
    assert np.abs(packet.energy - init_energy) > 1e-7


@pytest.mark.parametrize(
    "line_interaction_type",
    [
        LineInteractionType.SCATTER,
        LineInteractionType.DOWNBRANCH,
        LineInteractionType.MACROATOM,
    ],
)
def test_nonhomologous_line_scatter(
    line_interaction_type,
    packet,
    verysimple_time_explosion,
    verysimple_opacity_state,
    verysimple_numba_active_radial_1d_geometry,
):
    """
    Analogous to
    tardis/transport/montecarlo/tests/test_interaction.py@test_line_scatter
    """
    init_mu = packet.mu
    init_nu = packet.nu
    init_energy = packet.energy
    full_relativity = False
    packet.initialize_line_id(
        verysimple_opacity_state, verysimple_time_explosion, full_relativity
    )

    nonhomologous_line_scatter_event(
        packet,
        verysimple_numba_active_radial_1d_geometry,
        line_interaction_type,
        verysimple_opacity_state,
        enable_full_relativity=False,
    )

    assert np.abs(packet.mu - init_mu) > 1e-7
    assert np.abs(packet.nu - init_nu) > 1e-7
    assert np.abs(packet.energy - init_energy) > 1e-7


@pytest.mark.parametrize(
    ["test_packet", "expected"],
    [
        (
            {
                "mu": 0.8599443103322428,
                "emission_line_id": 1000,
                "energy": 0.9114437898710559,
                "nu": 0.0,
            },
            {"mu": 0.8599443103322428, "energy": 0.9114437898710559},
        ),
        (
            {
                "mu": -0.6975116557422458,
                "emission_line_id": 2000,
                "energy": 0.8803098648913266,
            },
            {"mu": -0.6975116557422458, "energy": 0.8803098648913266},
        ),
        (
            {
                "mu": -0.7115661419975774,
                "emission_line_id": 0,
                "energy": 0.8800385929341252,
            },
            {"mu": -0.7115661419975774, "energy": 0.8800385929341252},
        ),
    ],
)
def test_nonhomologous_line_emission(
    packet,
    verysimple_time_explosion,
    verysimple_opacity_state,
    verysimple_numba_active_radial_1d_geometry,
    test_packet,
    expected,
):
    """
    Analogous to
    tardis/transport/montecarlo/tests/test_interaction.py@test_line_emission
    """
    emission_line_id = test_packet["emission_line_id"]
    packet.mu = test_packet["mu"]
    packet.energy = test_packet["energy"]
    full_relativity = False
    packet.initialize_line_id(
        verysimple_opacity_state, verysimple_time_explosion, full_relativity
    )

    nonhomologous_line_emission(
        packet,
        emission_line_id,
        verysimple_numba_active_radial_1d_geometry,
        verysimple_opacity_state,
        full_relativity,
    )

    assert packet.next_line_id == emission_line_id + 1
    npt.assert_almost_equal(packet.mu, expected["mu"])
    npt.assert_almost_equal(packet.energy, expected["energy"])


def test_nonhomologous_calculate_sobolev_line_opacity(
    nb_simulation_verysimple,
    verysimple_numba_active_radial_1d_geometry,
    regression_data,
):
    """
    Analogous to
    tardis/opacities/tests/test_tau_sobolev.py@test_calculate_sobolev_line_opacity
    """
    legacy_plasma = nb_simulation_verysimple.plasma
    velocity_gradient = (
        verysimple_numba_active_radial_1d_geometry.velocity_gradient
        * u.Unit("1/s")
    )

    actual = nonhomologous_calculate_sobolev_line_opacity(
        legacy_plasma.atomic_data.lines,
        legacy_plasma.level_number_density,
        velocity_gradient,
        legacy_plasma.stimulated_emission_factor,
    )
    expected = regression_data.sync_dataframe(actual)
    pdt.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    "line_interaction_type,disable_line_scattering",
    [
        ("scatter", False),
        ("macroatom", False),
        ("macroatom", True),
        ("downbranch", False),
        ("downbranch", True),
    ],
)
def test_nonhomologous_opacity_solver(
    nb_simulation_verysimple,
    verysimple_numba_active_radial_1d_geometry,
    line_interaction_type,
    disable_line_scattering,
):
    """
    Analogous to
    tardis/opacities/tests/test_opacity_solver.py@test_opacity_solver
    """
    legacy_plasma = nb_simulation_verysimple.plasma
    velocity_gradient = (
        verysimple_numba_active_radial_1d_geometry.velocity_gradient
        * u.Unit("1/s")
    )

    solver = NonhomologousOpacitySolver(
        velocity_gradient=velocity_gradient,
        line_interaction_type=line_interaction_type,
        disable_line_scattering=disable_line_scattering,
    )
    actual = solver.legacy_solve(legacy_plasma)

    pdt.assert_series_equal(
        actual.electron_density, legacy_plasma.electron_densities
    )
    pdt.assert_series_equal(
        actual.line_list_nu, legacy_plasma.atomic_data.lines.nu
    )
    if not disable_line_scattering:
        pdt.assert_frame_equal(actual.tau_sobolev, legacy_plasma.tau_sobolevs)
    if line_interaction_type == "scatter":
        pass
    else:
        macro_atom_state = LegacyMacroAtomSolver().solve(
            legacy_plasma.j_blues,
            legacy_plasma.atomic_data,
            actual.tau_sobolev,
            legacy_plasma.stimulated_emission_factor,
            beta_sobolev=actual.beta_sobolev,
        )
        pdt.assert_frame_equal(
            macro_atom_state.transition_probabilities,
            legacy_plasma.transition_probabilities,
        )
        npt.assert_allclose(
            macro_atom_state.line2macro_level_upper,
            legacy_plasma.atomic_data.lines_upper2macro_reference_idx,
        )
        pdt.assert_series_equal(
            macro_atom_state.macro_block_edge_index,
            legacy_plasma.atomic_data.macro_atom_references["block_references"],
        )
        pdt.assert_series_equal(
            macro_atom_state.transition_type,
            legacy_plasma.atomic_data.macro_atom_data["transition_type"],
        )
        pdt.assert_series_equal(
            macro_atom_state.destination_level_id,
            legacy_plasma.atomic_data.macro_atom_data["destination_level_idx"],
        )
        pdt.assert_series_equal(
            macro_atom_state.transition_line_id,
            legacy_plasma.atomic_data.macro_atom_data["lines_idx"],
        )
