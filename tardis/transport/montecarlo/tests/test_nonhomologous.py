from collections.abc import Callable
from decimal import Decimal, localcontext
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from numpy.typing import NDArray

from tardis import constants as const
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.model.geometry.radial1d_homologous import (
    NumbaHomologousRadial1DGeometry,
)
from tardis.opacities.macro_atom.macroatom_solver import LegacyMacroAtomSolver
from tardis.opacities.opacity_state import OpacityState
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
from tardis.transport.montecarlo.configuration.base import (
    MonteCarloConfiguration,
)
from tardis.transport.montecarlo.configuration.constants import C_SPEED_OF_LIGHT
from tardis.transport.montecarlo.estimators.estimators_bulk import (
    EstimatorsBulk,
    init_estimators_bulk,
)
from tardis.transport.montecarlo.estimators.estimators_line import (
    EstimatorsLine,
    init_estimators_line,
)
from tardis.transport.montecarlo.modes.classic.packet_propagation import (
    packet_propagation as propagate_homologous,
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
from tardis.transport.montecarlo.modes.nonhomologous.packet_propagation import (
    packet_propagation as propagate_nonhomologous,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    trace_packet as nonhomologous_trace_packet,
)
from tardis.transport.montecarlo.modes.nonhomologous.rad_packet_transport import (
    update_line_ids,
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
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)


@pytest.fixture
def homologous_limit_case(
    homologous_radial_1d_geometry: NumbaHomologousRadial1DGeometry,
    classic_opacity_state: OpacityStateNumba,
) -> SimpleNamespace:
    """Provide an isolated line in a linear v=r/t field in two shells."""
    homologous = homologous_radial_1d_geometry
    nonhomologous = NumbaRadial1DGeometry(
        homologous.r_inner,
        homologous.r_outer,
        homologous.v_inner,
        homologous.v_outer,
    )
    classic_opacity_state.electron_density[:] = 0.0
    classic_opacity_state.tau_sobolev[:] = 0.0
    classic_opacity_state.tau_sobolev[0, :] = 100.0
    classic_opacity_state.line_list_nu[:] = [3.9999e14, 1.0e9]
    radius, mu, frequency = 7.5e14, 0.3, 4.0e14
    ct = const.c.cgs.value * homologous.time_explosion
    lab_frequency = frequency / (1.0 - radius * mu / ct)
    # Decimal solves the unsquared linear equation, not the production quartic.
    with localcontext() as context:
        context.prec = 80
        analytic_distance_decimal = Decimal.from_float(ct) * (
            1
            - Decimal.from_float(3.9999e14) / Decimal.from_float(lab_frequency)
        ) - Decimal.from_float(radius) * Decimal.from_float(mu)
    analytic_distance = float(analytic_distance_decimal)
    neighbourhood = 0.1 * (homologous.r_outer[0] - homologous.r_inner[0])
    maximum_projected_radius = (
        abs(radius * mu) + analytic_distance + neighbourhood
    )
    maximum_projected_beta = maximum_projected_radius / ct
    # IEEE binary64: three operations in (r*mu+s)/ct give gamma_3*A.
    # The two values near unity each round by <=2^-54; their final
    # subtraction is exact by Sterbenz's lemma (both lie between 1/2 and 1).
    unit_roundoff = np.finfo(np.float64).eps / 2.0
    gamma_three = 3.0 * unit_roundoff / (1.0 - 3.0 * unit_roundoff)
    assert 0.5 < 1.0 - maximum_projected_beta < 1.0
    assert 0.5 < 3.9999e14 / lab_frequency < 1.0
    residual_bound = gamma_three * maximum_projected_beta + unit_roundoff
    minimum_derivative = (
        1.0 / ct
    )  # |F'(s)|, in cm^-1 throughout the neighbourhood.
    distance_tolerance = 10.0 * residual_bound / minimum_derivative
    return SimpleNamespace(
        homologous=homologous,
        nonhomologous=nonhomologous,
        opacity=classic_opacity_state,
        time=homologous.time_explosion,
        radius=radius,
        mu=mu,
        lab_frequency=lab_frequency,
        ct=ct,
        analytic_distance=analytic_distance,
        analytic_distance_decimal=analytic_distance_decimal,
        neighbourhood=neighbourhood,
        residual_bound=residual_bound,
        distance_tolerance=distance_tolerance,
    )


def test_nonhomologous_resonance_matches_homologous_limit(
    homologous_limit_case: SimpleNamespace,
    line_estimators: EstimatorsLine,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """Verify the root and cursor against the homologous solution (verification).

    Claim: v=r/t, dv/dr=1/t gives the analytic resonance and ordered cursor.
    Regime: oblique outward ray, isolated opaque line well inside shell zero,
    with a red sentinel well beyond the outer boundary. Verification: solve
    F(s)=1-(r*mu+s)/(ct)-nu_line/nu_lab=0, independently of the quartic.
    On the neighbourhood s_exact +/- 0.1 shell width, |F'|=1/(ct) cm^-1.
    The independent 80-digit Decimal root uses the binary64 physical inputs.
    Over this neighbourhood, binary64 evaluation of F has absolute error
    <=gamma_3*A+u, A=max(|r*mu|+s)/(ct), u=2^-53: three operations form A,
    each near-unity subtraction/division rounds by <=2^-54, and the final
    subtraction is exact by Sterbenz. The resulting residual bound is about
    1.111e-16; ten times this bound / min|F'| gives about 1732 cm, selected
    without inspecting the production root. Exact event/cursor/shell assertions
    detect wrong-root selection, a wrong gradient or an off-by-one cursor.
    """
    case = homologous_limit_case
    radius, mu, lab_frequency = case.radius, case.mu, case.lab_frequency
    ct = case.ct
    packet = RPacket(radius, mu, lab_frequency, 0.9, 42, 0)
    packet.next_line_id, packet.prev_line_id = 0, -1
    line_frequency = case.opacity.line_list_nu[0]
    analytic_distance = case.analytic_distance
    distance_tolerance = case.distance_tolerance
    neighbourhood = case.neighbourhood
    boundary_distance = (
        np.sqrt(case.homologous.r_outer[0] ** 2 - radius**2 * (1.0 - mu**2))
        - radius * mu
    )
    assert 0.0 < analytic_distance - neighbourhood
    assert analytic_distance + neighbourhood < boundary_distance
    actual_distance = calculate_distance_line_nonhomologous(
        packet, case.nonhomologous, line_frequency
    )
    npt.assert_allclose(
        actual_distance, analytic_distance, rtol=0, atol=distance_tolerance
    )
    # Also check the unsquared physical equation, which rejects extraneous roots.
    with localcontext() as context:
        context.prec = 80
        residual = (
            case.analytic_distance_decimal - Decimal.from_float(actual_distance)
        ) / Decimal.from_float(ct)
    assert abs(residual) <= Decimal.from_float(10.0 * case.residual_bound)
    npt.assert_allclose(
        case.nonhomologous.velocity_gradient,
        1.0 / case.time,
        rtol=1.0e-12,
        atol=0,
    )
    set_seed_fixture(42)
    distance, event, delta_shell = nonhomologous_trace_packet(
        packet,
        case.nonhomologous,
        case.opacity,
        line_estimators,
        0.0,
        False,
        False,
    )
    assert event == InteractionType.LINE
    assert packet.current_shell_id == 0
    assert packet.next_line_id == 0
    assert packet.prev_line_id == -1
    assert delta_shell == 1
    npt.assert_allclose(
        distance, analytic_distance, rtol=0, atol=distance_tolerance
    )
    expected_energy = packet.energy * line_frequency / lab_frequency
    # A root displacement ds changes E_comov by E_lab*ds/(ct), not by ds erg.
    energy_tolerance = packet.energy * distance_tolerance / ct
    npt.assert_allclose(
        line_estimators.energy_deposition_line_rate[0, 0],
        expected_energy,
        rtol=1.0e-12,
        atol=energy_tolerance,
    )
    npt.assert_allclose(
        line_estimators.mean_intensity_blueward[0, 0],
        expected_energy / lab_frequency,
        rtol=1.0e-12,
        atol=energy_tolerance / lab_frequency,
    )


def test_nonhomologous_packet_matches_homologous_limit(
    homologous_limit_case: SimpleNamespace,
    montecarlo_configuration: MonteCarloConfiguration,
    set_seed_fixture: Callable[[int], None],
) -> None:
    """Characterize seeded end-to-end agreement in the homologous limit.

    Claim: the current propagators agree on final packet, event/cursor history
    and bulk/line estimators for v=r/t and dv/dr=1/t. Regime: one isolated
    opaque line, radiative scattering, electron-free two-shell ejecta, seed 42.
    Verification: metamorphic comparison with homologous transport; the root
    is independently verified in the companion test. Evidence: characterization
    (not an independent validation of either entire transport solver).
    The companion fixture derives delta~=1732 cm from an independent IEEE
    residual bound and an 80-digit unsquared root. Two roots can differ by
    Delta=2*delta. For the fixed outward emission cosine z>0, |dr_line/ds|<=1,
    |dD_in/ds|=1/(ct), |dD_out/dr|=z/(ct). Thus packet E proportional to
    D_in/D_out and nu proportional to 1/D_out have separately bounded relative
    errors a+b+ab and b. At a boundary R, |dmu_R/dr_line| is bounded by
    r_max*(1-z^2)/(R^2*z). Segment length derivatives are 1, <=1/z, and
    <=r_max*(1-z^2)/z*(1/R_inner+1/R_outer) for the three positive path lengths.
    Each bulk E*D*length or E*nu*D^2*length product uses these relative bounds;
    a shell sum uses its largest segment bound. No distance in cm is reused
    as an energy/frequency tolerance. Floating point allowance is 1e-12.
    Exact discrete outcomes detect wrong cursor/process/dispatcher wiring.
    """
    case = homologous_limit_case
    packets, trackers, bulks, lines = [], [], [], []
    for geometry, propagate in [
        (case.homologous, propagate_homologous),
        (case.nonhomologous, propagate_nonhomologous),
    ]:
        packet = RPacket(7.5e14, 0.3, 4.0e14, 0.9, 42, 0)
        tracker = TrackerFull(16, 2)
        bulk, line = init_estimators_bulk(2), init_estimators_line((2, 2))
        set_seed_fixture(42)
        propagate(
            packet,
            geometry,
            case.time,
            case.opacity,
            bulk,
            line,
            tracker,
            montecarlo_configuration,
        )
        tracker.finalize()
        packets.append(packet)
        trackers.append(tracker)
        bulks.append(bulk)
        lines.append(line)
    classic, nonhomologous = packets
    ct = case.ct
    distance_difference = 2.0 * case.distance_tolerance
    # The second MT19937 draw is the isotropic emission direction, common
    # to both paths; the first draw selects the opaque line interaction.
    emission_mu = 2.0 * np.random.RandomState(42).random_sample(2)[1] - 1.0
    assert emission_mu > 0.0
    line_radius = np.sqrt(
        case.radius**2
        + case.analytic_distance**2
        + 2.0 * case.radius * case.analytic_distance * case.mu
    )
    maximum_line_radius = line_radius + case.neighbourhood
    inner_boundary, outer_boundary = case.homologous.r_outer
    assert maximum_line_radius < inner_boundary
    minimum_doppler = 1.0 - inner_boundary / ct
    incoming_doppler_fraction = distance_difference / (ct * minimum_doppler)
    outgoing_doppler_fraction = emission_mu * incoming_doppler_fraction
    energy_fraction = (
        incoming_doppler_fraction
        + outgoing_doppler_fraction
        + incoming_doppler_fraction * outgoing_doppler_fraction
    )
    frequency_fraction = outgoing_doppler_fraction
    final_mu_error = (
        maximum_line_radius
        * (1.0 - emission_mu**2)
        * distance_difference
        / (outer_boundary**2 * emission_mu)
    )
    for packet in packets:
        assert packet.status == PacketStatus.EMITTED
        assert packet.current_shell_id == 1
        assert packet.next_line_id == 1
        npt.assert_allclose(
            packet.r, case.homologous.r_outer[-1], rtol=1.0e-12, atol=0
        )
    npt.assert_allclose(
        nonhomologous.mu, classic.mu, rtol=0, atol=final_mu_error + 1.0e-12
    )
    npt.assert_allclose(
        nonhomologous.energy,
        classic.energy,
        rtol=energy_fraction + 1.0e-12,
        atol=0,
    )
    npt.assert_allclose(
        nonhomologous.nu, classic.nu, rtol=frequency_fraction + 1.0e-12, atol=0
    )
    for attribute in [
        "interaction_type",
        "before_shell_id",
        "after_shell_id",
        "line_absorb_id",
        "line_emit_id",
        "next_line_id",
        "status",
    ]:
        npt.assert_array_equal(
            getattr(trackers[1], attribute), getattr(trackers[0], attribute)
        )
    npt.assert_array_equal(
        trackers[0].interaction_type,
        [
            InteractionType.BOUNDARY,
            InteractionType.LINE,
            InteractionType.BOUNDARY,
            InteractionType.BOUNDARY,
            InteractionType.BOUNDARY,
        ],
    )
    npt.assert_allclose(
        trackers[1].radius,
        trackers[0].radius,
        rtol=1.0e-12,
        atol=distance_difference,
    )
    for tracker in trackers:
        npt.assert_allclose(
            tracker.radius[1],
            line_radius,
            rtol=1.0e-12,
            atol=case.distance_tolerance,
        )
        npt.assert_allclose(
            tracker.after_mu[1], emission_mu, rtol=0, atol=1.0e-15
        )
    npt.assert_allclose(
        lines[1].energy_deposition_line_rate,
        lines[0].energy_deposition_line_rate,
        rtol=incoming_doppler_fraction + 1.0e-12,
        atol=0,
    )
    npt.assert_allclose(
        lines[1].mean_intensity_blueward,
        lines[0].mean_intensity_blueward,
        rtol=incoming_doppler_fraction + 1.0e-12,
        atol=0,
    )
    # Three segments: source->line, line->shell edge, shell edge->outer edge.
    impact_squared = line_radius**2 * (1.0 - emission_mu**2)
    edge_projections = np.sqrt(case.homologous.r_outer**2 - impact_squared)
    segment_lengths = np.array(
        [
            case.analytic_distance,
            edge_projections[0] - line_radius * emission_mu,
            edge_projections[1] - edge_projections[0],
        ]
    )
    length_derivatives = np.array(
        [
            1.0,
            1.0 / emission_mu,
            maximum_line_radius
            * (1.0 - emission_mu**2)
            / emission_mu
            * (1.0 / inner_boundary + 1.0 / outer_boundary),
        ]
    )
    minimum_lengths = (
        segment_lengths - length_derivatives * case.distance_tolerance
    )
    assert np.all(minimum_lengths > 0.0)
    length_fractions = (
        length_derivatives * distance_difference / minimum_lengths
    )
    doppler_derivatives = np.array(
        [
            0.0,
            emission_mu / ct,
            maximum_line_radius
            * (1.0 - emission_mu**2)
            / (ct * inner_boundary * emission_mu),
        ]
    )
    doppler_fractions = (
        doppler_derivatives * distance_difference / minimum_doppler
    )
    energy_fractions = np.array([0.0, energy_fraction, energy_fraction])
    frequency_fractions = np.array(
        [0.0, frequency_fraction, frequency_fraction]
    )
    # Positive product perturbations: prod(1+relative_error)-1, evaluated
    # with log1p/expm1 to avoid cancellation at ~1e-11 relative errors.
    intensity_fractions = np.expm1(
        np.log1p(energy_fractions)
        + np.log1p(doppler_fractions)
        + np.log1p(length_fractions)
    )
    mean_frequency_fractions = np.expm1(
        np.log1p(energy_fractions)
        + np.log1p(frequency_fractions)
        + 2.0 * np.log1p(doppler_fractions)
        + np.log1p(length_fractions)
    )
    for shell_id, segment_ids in [(0, [0, 1]), (1, [2])]:
        npt.assert_allclose(
            bulks[1].mean_intensity_total[shell_id],
            bulks[0].mean_intensity_total[shell_id],
            rtol=np.max(intensity_fractions[segment_ids]) + 1.0e-12,
            atol=0,
        )
        npt.assert_allclose(
            bulks[1].mean_frequency[shell_id],
            bulks[0].mean_frequency[shell_id],
            rtol=np.max(mean_frequency_fractions[segment_ids]) + 1.0e-12,
            atol=0,
        )


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


@pytest.fixture
def distance_solver_geometry(line_search_packet: RPacket) -> NumbaRadial1DGeometry:
    """Return a narrow shell with a 2.5e5 cm/s velocity difference."""
    inner_radius = line_search_packet.r
    inner_velocity = inner_radius * 1.0e-7
    return NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([inner_radius + 1.0e12]),
        np.asarray([inner_velocity]),
        np.asarray([inner_velocity + 2.5e5]),
    )


@pytest.fixture
def gradient_reversal_setup(
    classic_opacity_state: OpacityStateNumba,
    line_search_packet: RPacket,
) -> tuple[NumbaRadial1DGeometry, NDArray[np.float64], NDArray[np.float64], float]:
    """Arrange a line with two analytic roots on opposite sides of a reversal.

    The unsquared projected velocity at x=0.25b and x=1.25b sets the linear
    field and line frequency. Its derivative supplies independent gradients.
    """
    impact_parameter = 1.0e14  # cm
    first_x, second_x = 0.25 * impact_parameter, 1.25 * impact_parameter
    velocity_intercept = 2.0e7  # cm / s
    first_radius = np.hypot(first_x, impact_parameter)
    second_radius = np.hypot(second_x, impact_parameter)
    velocity_gradient = velocity_intercept * (
        second_x / second_radius - first_x / first_radius
    ) / (first_x - second_x)
    inner_radius, outer_radius = 0.9 * impact_parameter, 1.7 * impact_parameter
    geometry = NumbaRadial1DGeometry(
        np.asarray([inner_radius]),
        np.asarray([outer_radius]),
        np.asarray([velocity_gradient * inner_radius + velocity_intercept]),
        np.asarray([velocity_gradient * outer_radius + velocity_intercept]),
    )
    first_projected_velocity = (
        velocity_gradient * first_x
        + velocity_intercept * first_x / first_radius
    )
    line_frequency = line_search_packet.nu * (
        1.0 - first_projected_velocity / C_SPEED_OF_LIGHT
    )
    classic_opacity_state.line_list_nu = np.asarray([line_frequency, 1.0e14])
    turning_radius = (
        -velocity_intercept * impact_parameter**2 / velocity_gradient
    ) ** (1.0 / 3.0)
    expected_zero = np.sqrt(turning_radius**2 - impact_parameter**2)
    first_gradient = velocity_gradient + (
        velocity_intercept * impact_parameter**2 / first_radius**3
    )
    second_gradient = velocity_gradient + (
        velocity_intercept * impact_parameter**2 / second_radius**3
    )
    return (
        geometry,
        np.asarray([first_x, second_x]),
        np.asarray([first_gradient, second_gradient]),
        expected_zero,
    )


@pytest.fixture
def close_line_setup(
    classic_opacity_state: OpacityStateNumba,
    line_search_packet: RPacket,
    line_search_velocity_parameters: tuple[float, float],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Arrange two rest frequencies using the unsquared ray equation."""
    velocity_gradient, velocity_intercept = line_search_velocity_parameters
    packet = line_search_packet
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
    classic_opacity_state.line_list_nu = line_frequencies
    return np.asarray([first_distance, second_distance]), line_frequencies


@pytest.fixture
def frequency_changing_emission_setup(
    classic_opacity_state: OpacityStateNumba,
    line_search_packet: RPacket,
    line_search_velocity_parameters: tuple[float, float],
) -> tuple[float, float, float]:
    """Arrange absorption at 1e12 cm and a repeated crossing after re-emission.

    The unsquared ray equation fixes the incoming lab frequency and a higher
    emitted line frequency whose next crossing is 3e12 cm after emission.
    """
    velocity_gradient, velocity_intercept = line_search_velocity_parameters
    packet = line_search_packet
    initial_radius, initial_mu = packet.r, packet.mu
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
    packet.nu = incoming_lab_frequency
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
    classic_opacity_state.line_list_nu = np.asarray(
        [emitted_frequency, absorbed_frequency]
    )
    classic_opacity_state.tau_sobolev[1, 0] = 1.0e6
    return first_line_distance, repeated_line_distance, initial_comov_frequency


@pytest.mark.parametrize(
    "line_search_packet", [{"r": 1.0e14, "mu": 0.0}], indirect=True
)
def test_nonhomologous_distance_preserves_near_line_center_root(
    line_search_packet: RPacket,
    distance_solver_geometry: NumbaRadial1DGeometry,
) -> None:
    """Claim: The distance solver selects the physical near-line-center root.

    Regime: A tangential ray in a narrow linear-velocity shell, where squaring
    the resonance equation produces a nearby nonphysical root.
    Verification: The distance matches the reference physical root
    0.14167578672128847 times the 1e12 cm shell width.
    """
    shell_width = 1.0e12
    geometry, packet = distance_solver_geometry, line_search_packet
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


@pytest.mark.parametrize(
    "line_search_packet", [{"r": 2.0e14, "mu": 0.0}], indirect=True
)
def test_nonhomologous_distance_rejects_sign_reversed_root(
    line_search_packet: RPacket,
    distance_solver_geometry: NumbaRadial1DGeometry,
) -> None:
    """Claim: Squaring the resonance equation must not select a spurious root.

    Regime: A tangential ray in a linear-velocity shell with a closer,
    sign-reversed quartic root.
    Verification: Construct the line frequency from the unsquared equation
    at a prescribed projected position and recover that position as distance.
    """
    shell_width = 1.0e12
    shell_velocity_difference = 2.5e5
    geometry, packet = distance_solver_geometry, line_search_packet

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


@pytest.mark.parametrize(
    "line_search_packet", [{"r": 1.0e14, "mu": 0.0}], indirect=True
)
def test_nonhomologous_line_reappears_after_gradient_reversal(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
    line_search_packet: RPacket,
    gradient_reversal_setup: tuple[
        NumbaRadial1DGeometry, NDArray[np.float64], NDArray[np.float64], float
    ],
) -> None:
    """Claim: A line contributes again after an interior gradient reversal.

    Regime: A nonradial ray in a linear-velocity shell with two resonances
    and negligible electron opacity; line scattering is disabled.
    Verification: The unsquared ray equation fixes both resonance positions
    and gradients independently of the solver. Compare the turning point,
    both distances, and the summed energy and blueward estimator contributions.
    """
    geometry, positions, gradients, expected_zero = gradient_reversal_setup
    packet = line_search_packet
    first_x, second_x = positions
    first_gradient, second_gradient = gradients
    line_frequency = classic_opacity_state.line_list_nu[0]

    boundary_distance, _ = calculate_distance_boundary(
        packet.r, packet.mu, geometry.r_inner[0], geometry.r_outer[0]
    )
    first_zero, _, zero_count = calculate_projected_gradient_zero_distances(
        packet, geometry, boundary_distance
    )
    assert zero_count == 1
    npt.assert_allclose(first_zero, expected_zero, rtol=0.0, atol=100.0)

    first_distance = calculate_distance_line_nonhomologous(
        packet, geometry, line_frequency, 0.0, first_zero
    )
    second_distance = calculate_distance_line_nonhomologous(
        packet, geometry, line_frequency, first_zero, boundary_distance
    )
    # Adjacent representable frequencies near 1e15 Hz correspond to less than
    # 100 cm in resonance position at these gradients.
    npt.assert_allclose(
        [first_distance, second_distance],
        [first_x, second_x],
        rtol=0.0,
        atol=100.0,
    )

    update_line_ids(
        packet,
        classic_opacity_state.line_list_nu,
        calculate_comoving_frequency_nonhomologous(packet, geometry, 0.0),
    )
    _, interaction_type, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        True,
    )
    assert interaction_type == InteractionType.BOUNDARY

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


def test_nonhomologous_search_resolves_two_close_lines(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
    line_search_geometry: NumbaRadial1DGeometry,
    line_search_packet: RPacket,
    close_line_setup: tuple[NDArray[np.float64], NDArray[np.float64]],
) -> None:
    """Claim: Two nearby distinct lines each contribute on a monotonic ray.

    Regime: A nonradial ray in a positive-gradient shell, with resonances
    separated by 1e7 cm, negligible electron opacity, and no line scattering.
    Verification: Frequencies derived from the unsquared ray equation must
    recover the prescribed distances and separate Doppler-weighted energies.
    """
    geometry, packet = line_search_geometry, line_search_packet
    positions, line_frequencies = close_line_setup
    first_distance, second_distance = positions
    assert line_frequencies[0] > line_frequencies[1]

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

    update_line_ids(
        packet,
        classic_opacity_state.line_list_nu,
        calculate_comoving_frequency_nonhomologous(packet, geometry, 0.0),
    )
    _, interaction_type, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        True,
    )
    assert interaction_type == InteractionType.BOUNDARY
    npt.assert_allclose(
        line_estimators.energy_deposition_line_rate[:, 0],
        line_frequencies / packet.nu,
        rtol=1.0e-10,
    )


def test_nonhomologous_frequency_changing_emission_starts_new_line_history(
    classic_opacity_state: OpacityStateNumba,
    line_estimators: EstimatorsLine,
    bulk_estimators: EstimatorsBulk,
    line_search_geometry: NumbaRadial1DGeometry,
    line_search_packet: RPacket,
    frequency_changing_emission_setup: tuple[float, float, float],
) -> None:
    """Claim: Emission at a higher rest frequency starts a new line history.

    Regime: Two lines in a monotonic nonhomologous flow, with absorption in
    the lower-frequency line followed by emission in the higher-frequency one.
    Verification: The unsquared ray equation prescribes crossings at 1e12 cm
    and 3e12 cm after emission. Check both events, the emitted cursor pair,
    and the repeated line's new Doppler-weighted energy contribution.
    """
    geometry, packet = line_search_geometry, line_search_packet
    first_line_distance, repeated_line_distance, initial_comov_frequency = (
        frequency_changing_emission_setup
    )
    emitted_frequency, absorbed_frequency = classic_opacity_state.line_list_nu
    assert emitted_frequency > initial_comov_frequency > absorbed_frequency

    update_line_ids(
        packet, classic_opacity_state.line_list_nu, initial_comov_frequency
    )
    first_distance, first_event, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        False,
    )
    assert first_event == InteractionType.LINE
    assert packet.next_line_id == 1
    npt.assert_allclose(first_distance, first_line_distance, rtol=0.0, atol=100.0)
    first_deposition = line_estimators.energy_deposition_line_rate[1, 0]

    move_r_packet(packet, first_distance, geometry, bulk_estimators, False)
    nonhomologous_line_emission(
        packet,
        0,
        geometry,
        classic_opacity_state,
        False,
    )
    # Emission excludes its own ID and makes line 1 eligible on the new flight.
    assert (packet.next_line_id, packet.prev_line_id) == (1, -1)

    second_distance, second_event, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        classic_opacity_state,
        line_estimators,
        1.0e-30,
        False,
        False,
    )
    assert second_event == InteractionType.LINE
    assert packet.next_line_id == 1
    npt.assert_allclose(
        second_distance, repeated_line_distance, rtol=0.0, atol=100.0
    )
    npt.assert_allclose(
        line_estimators.energy_deposition_line_rate[1, 0] - first_deposition,
        packet.energy * absorbed_frequency / packet.nu,
        rtol=1.0e-10,
    )
    assert line_estimators.energy_deposition_line_rate[0, 0] == 0.0


@pytest.fixture(params=[1, -1], ids=["redshift", "blueshift"])
def blend_line_id_step(request: pytest.FixtureRequest) -> int:
    """Select the drift sign in a finite-gradient radial blend scenario."""
    return request.param


@pytest.fixture
def blend_velocity_parameters(blend_line_id_step: int) -> tuple[float, float]:
    """Return a linear positive-velocity field with abs(G)=8e-7 s^-1."""
    return (8.0e-7, 2.0e7) if blend_line_id_step == 1 else (-8.0e-7, 2.6e8)


@pytest.fixture
def blend_packet(line_search_packet: RPacket) -> RPacket:
    """Launch the shared unit-energy packet on a radial outgoing ray."""
    line_search_packet.mu = 1.0
    return line_search_packet


@pytest.fixture
def blend_geometry(
    blend_velocity_parameters: tuple[float, float],
) -> NumbaRadial1DGeometry:
    """Return one shell for the common blend velocity field."""
    gradient, intercept = blend_velocity_parameters
    inner = np.array([1.0e14])
    outer = np.array([2.0e14])
    return NumbaRadial1DGeometry(
        inner, outer, gradient * inner + intercept, gradient * outer + intercept
    )


@pytest.fixture
def blend_resonance_distance(request: pytest.FixtureRequest) -> float:
    """Prescribe an independent radial resonance position, in centimeters."""
    return getattr(request, "param", 0.0)


@pytest.fixture
def blend_opacity_state(
    classic_opacity_state: OpacityStateNumba,
    blend_packet: RPacket,
    blend_velocity_parameters: tuple[float, float],
    blend_resonance_distance: float,
) -> OpacityStateNumba:
    """Give three transition IDs the frequency of one analytic resonance."""
    gradient, intercept = blend_velocity_parameters
    resonance_velocity = (
        gradient * (blend_packet.r + blend_resonance_distance) + intercept
    )
    frequency = blend_packet.nu * (1.0 - resonance_velocity / C_SPEED_OF_LIGHT)
    classic_opacity_state.line_list_nu = np.full(3, frequency)
    classic_opacity_state.tau_sobolev = np.zeros((3, 2))
    classic_opacity_state.sobolev_optical_depth_coefficient = np.zeros((3, 2))
    return classic_opacity_state


@pytest.fixture
def blend_estimators() -> EstimatorsLine:
    """Keep a separate estimator channel for every blend member and shell."""
    return EstimatorsLine(np.zeros((3, 2)), np.zeros((3, 2)))


@pytest.fixture
def blend_trace(
    monkeypatch: pytest.MonkeyPatch,
    blend_packet: RPacket,
    blend_geometry: NumbaRadial1DGeometry,
    blend_opacity_state: OpacityStateNumba,
) -> Callable[..., tuple[float, InteractionType, int]]:
    """Initialize the real tracer's history and fix event optical depth to one."""
    update_line_ids(
        blend_packet,
        blend_opacity_state.line_list_nu,
        calculate_comoving_frequency_nonhomologous(
            blend_packet, blend_geometry, 0.0
        ),
    )
    monkeypatch.setattr(np.random, "random", lambda: np.exp(-1.0))
    return getattr(
        nonhomologous_trace_packet, "py_func", nonhomologous_trace_packet
    )


def test_nonhomologous_emission_processes_eligible_blend_at_zero_distance(
    blend_line_id_step: int,
    blend_velocity_parameters: tuple[float, float],
    blend_packet: RPacket,
    blend_geometry: NumbaRadial1DGeometry,
    blend_opacity_state: OpacityStateNumba,
    blend_estimators: EstimatorsLine,
    blend_trace: Callable[..., tuple[float, InteractionType, int]],
) -> None:
    """Claim: Emission excludes its own ID while leaving a blend twin eligible.

    Regime: A radial ray with finite redward or blueward drift, three
    equal-frequency transitions, and emission in the middle transition.
    Verification: The directional twin scatters at exactly zero distance and
    receives the analytic comoving energy and blueward contributions; all
    other transition channels remain zero.
    """
    gradient, _ = blend_velocity_parameters
    blend_opacity_state.sobolev_optical_depth_coefficient[:] = 2.0 * abs(
        gradient
    )
    nonhomologous_line_emission(
        blend_packet, 1, blend_geometry, blend_opacity_state, False
    )
    eligible_line_id = 1 + blend_line_id_step
    expected_energy = (
        blend_packet.energy
        * blend_opacity_state.line_list_nu[eligible_line_id]
        / blend_packet.nu
    )

    distance, event, _ = blend_trace(
        blend_packet,
        blend_geometry,
        blend_opacity_state,
        blend_estimators,
        1.0e-30,
        False,
        False,
    )

    assert event == InteractionType.LINE
    assert distance == 0.0
    assert blend_packet.next_line_id == eligible_line_id
    expected_deposition = np.zeros(3)
    expected_deposition[eligible_line_id] = expected_energy
    npt.assert_allclose(
        blend_estimators.energy_deposition_line_rate[:, 0],
        expected_deposition,
        rtol=1.0e-10,
        atol=0.0,
    )
    npt.assert_allclose(
        blend_estimators.mean_intensity_blueward[:, 0],
        expected_deposition / (blend_packet.nu * abs(gradient)),
        rtol=1.0e-10,
        atol=0.0,
    )


@pytest.mark.parametrize("blend_resonance_distance", [1.0e12], indirect=True)
def test_nonhomologous_blend_members_share_resonance_contributions(
    blend_velocity_parameters: tuple[float, float],
    blend_resonance_distance: float,
    blend_packet: RPacket,
    blend_geometry: NumbaRadial1DGeometry,
    blend_opacity_state: OpacityStateNumba,
    blend_estimators: EstimatorsLine,
    blend_trace: Callable[..., tuple[float, InteractionType, int]],
) -> None:
    """Claim: Every blend member contributes once at the common resonance.

    Regime: Three equal-frequency transitions on a radial ray in either drift
    direction, with finite gradient and line scattering disabled.
    Verification: The prescribed frequency gives a resonance at 1e12 cm.
    Each channel receives E_comov and E_comov / (nu_lab * abs(G)), and the
    trace reaches the boundary.
    """
    gradient, _ = blend_velocity_parameters
    line_frequency = blend_opacity_state.line_list_nu[0]
    expected_energy = blend_packet.energy * line_frequency / blend_packet.nu
    solved_distance = calculate_distance_line_nonhomologous(
        blend_packet, blend_geometry, line_frequency
    )
    npt.assert_allclose(
        solved_distance, blend_resonance_distance, rtol=0.0, atol=100.0
    )

    _, event, _ = blend_trace(
        blend_packet,
        blend_geometry,
        blend_opacity_state,
        blend_estimators,
        1.0e-30,
        False,
        True,
    )

    assert event == InteractionType.BOUNDARY
    npt.assert_allclose(
        blend_estimators.energy_deposition_line_rate[:, 0],
        np.full(3, expected_energy),
        rtol=1.0e-10,
    )
    npt.assert_allclose(
        blend_estimators.mean_intensity_blueward[:, 0],
        np.full(3, expected_energy / (blend_packet.nu * abs(gradient))),
        rtol=1.0e-10,
        atol=0.0,
    )


@pytest.mark.parametrize("blend_resonance_distance", [1.0e12], indirect=True)
def test_nonhomologous_blend_optical_depth_accumulates_in_order(
    blend_line_id_step: int,
    blend_velocity_parameters: tuple[float, float],
    blend_resonance_distance: float,
    blend_packet: RPacket,
    blend_geometry: NumbaRadial1DGeometry,
    blend_opacity_state: OpacityStateNumba,
    blend_estimators: EstimatorsLine,
    blend_trace: Callable[..., tuple[float, InteractionType, int]],
) -> None:
    """Claim: Blend optical depths accumulate in directional traversal order.

    Regime: A radial ray in either drift direction meets three transitions
    at 1e12 cm, with ordered depths 0.6, 0.6, and 2 and event depth one.
    Verification: Neither of the first two depths suffices alone; their sum
    selects the second member. Only the first two channels receive analytic
    energy and blueward contributions.
    """
    gradient, _ = blend_velocity_parameters
    # Neither of the first two depths can select a line without accumulation.
    depths = np.array([0.6, 0.6, 2.0])[::blend_line_id_step]
    blend_opacity_state.sobolev_optical_depth_coefficient[:, 0] = depths * abs(
        gradient
    )
    expected_energy = (
        blend_packet.energy
        * blend_opacity_state.line_list_nu[0]
        / blend_packet.nu
    )

    distance, event, _ = blend_trace(
        blend_packet,
        blend_geometry,
        blend_opacity_state,
        blend_estimators,
        1.0e-30,
        False,
        False,
    )

    assert event == InteractionType.LINE
    assert blend_packet.next_line_id == 1
    npt.assert_allclose(
        distance, blend_resonance_distance, rtol=0.0, atol=100.0
    )
    expected_deposition = np.full(3, expected_energy)
    expected_deposition[2 if blend_line_id_step == 1 else 0] = 0.0
    npt.assert_allclose(
        blend_estimators.energy_deposition_line_rate[:, 0],
        expected_deposition,
        rtol=1.0e-10,
        atol=0.0,
    )
    npt.assert_allclose(
        blend_estimators.mean_intensity_blueward[:, 0],
        expected_deposition / (blend_packet.nu * abs(gradient)),
        rtol=1.0e-10,
        atol=0.0,
    )


@pytest.mark.parametrize(
    "blend_resonance_distance", [2.5e13 - 200.0], indirect=True
)
def test_processed_line_survives_monotonic_boundary_return(
    blend_velocity_parameters: tuple[float, float],
    blend_packet: RPacket,
    blend_opacity_state: OpacityStateNumba,
    blend_estimators: EstimatorsLine,
    bulk_estimators: EstimatorsBulk,
    blend_trace: Callable[..., tuple[float, InteractionType, int]],
) -> None:
    """Claim: A continuous boundary return preserves processed-line history.

    Regime: A radial ray with unchanged drift sign crosses a blend 200 cm
    before a velocity-continuous shell edge, within the next trace's close
    band. Line scattering is disabled and electron opacity is negligible.
    Verification: Each transition receives one analytic energy and blueward
    contribution upstream and exactly zero contribution downstream.
    """
    gradient, intercept = blend_velocity_parameters
    inner = np.array([1.0e14, 1.5e14])
    outer = np.array([1.5e14, 2.0e14])
    geometry = NumbaRadial1DGeometry(
        inner, outer, gradient * inner + intercept, gradient * outer + intercept
    )
    expected_energy = (
        blend_packet.energy
        * blend_opacity_state.line_list_nu[0]
        / blend_packet.nu
    )
    distance, event, shell_step = blend_trace(
        blend_packet,
        geometry,
        blend_opacity_state,
        blend_estimators,
        1.0e-30,
        False,
        True,
    )
    assert event == InteractionType.BOUNDARY
    move_r_packet(blend_packet, distance, geometry, bulk_estimators, False)
    move_packet_across_shell_boundary(blend_packet, shell_step, 2)

    _, second_event, _ = blend_trace(
        blend_packet,
        geometry,
        blend_opacity_state,
        blend_estimators,
        1.0e-30,
        False,
        True,
    )

    assert second_event == InteractionType.BOUNDARY
    npt.assert_allclose(
        blend_estimators.energy_deposition_line_rate[:, 0],
        np.full(3, expected_energy),
        rtol=1.0e-10,
    )
    npt.assert_allclose(
        blend_estimators.mean_intensity_blueward[:, 0],
        np.full(3, expected_energy / (blend_packet.nu * abs(gradient))),
        rtol=1.0e-10,
        atol=0.0,
    )
    npt.assert_array_equal(
        blend_estimators.energy_deposition_line_rate[:, 1], 0.0
    )
    npt.assert_array_equal(blend_estimators.mean_intensity_blueward[:, 1], 0.0)


@pytest.mark.parametrize("blend_line_id_step", [1], indirect=True)
@pytest.mark.parametrize("blend_resonance_distance", [1.0e12], indirect=True)
def test_equal_frequency_different_line_emission_restarts_blend_search(
    blend_velocity_parameters: tuple[float, float],
    blend_packet: RPacket,
    blend_geometry: NumbaRadial1DGeometry,
    blend_opacity_state: OpacityStateNumba,
    blend_estimators: EstimatorsLine,
    bulk_estimators: EstimatorsBulk,
    blend_trace: Callable[..., tuple[float, InteractionType, int]],
) -> None:
    """Claim: Different-ID emission resets history even at equal rest frequency.

    Regime: A radial redshifting ray absorbs in ID 1 and emits in
    equal-frequency ID 0. ID 1 has depth two and each trace has event depth one.
    Verification: Both flights select ID 1; the second does so at zero
    distance. Compare each flight's analytic energy and blueward contributions
    and require no new contribution from emitting ID 0 on the second flight.
    """
    gradient, _ = blend_velocity_parameters
    blend_opacity_state.sobolev_optical_depth_coefficient[1, 0] = 2.0 * abs(
        gradient
    )
    expected_energy = (
        blend_packet.energy
        * blend_opacity_state.line_list_nu[1]
        / blend_packet.nu
    )
    distance, event, _ = blend_trace(
        blend_packet,
        blend_geometry,
        blend_opacity_state,
        blend_estimators,
        1.0e-30,
        False,
        False,
    )
    assert event == InteractionType.LINE
    assert blend_packet.next_line_id == 1
    expected_first_deposition = np.array([expected_energy, expected_energy, 0.0])
    npt.assert_allclose(
        blend_estimators.energy_deposition_line_rate[:, 0],
        expected_first_deposition,
        rtol=1.0e-10,
        atol=0.0,
    )
    npt.assert_allclose(
        blend_estimators.mean_intensity_blueward[:, 0],
        expected_first_deposition / (blend_packet.nu * abs(gradient)),
        rtol=1.0e-10,
        atol=0.0,
    )
    move_r_packet(
        blend_packet, distance, blend_geometry, bulk_estimators, False
    )
    nonhomologous_line_emission(
        blend_packet, 0, blend_geometry, blend_opacity_state, False
    )
    before_deposition = blend_estimators.energy_deposition_line_rate[
        :, 0
    ].copy()
    before_blueward = blend_estimators.mean_intensity_blueward[:, 0].copy()
    expected_energy = (
        blend_packet.energy
        * blend_opacity_state.line_list_nu[1]
        / blend_packet.nu
    )

    distance, event, _ = blend_trace(
        blend_packet,
        blend_geometry,
        blend_opacity_state,
        blend_estimators,
        1.0e-30,
        False,
        False,
    )

    assert event == InteractionType.LINE
    assert distance == 0.0
    assert blend_packet.next_line_id == 1
    expected_deposition = np.array([0.0, expected_energy, 0.0])
    npt.assert_allclose(
        blend_estimators.energy_deposition_line_rate[:, 0] - before_deposition,
        expected_deposition,
        rtol=1.0e-10,
        atol=0.0,
    )
    npt.assert_allclose(
        blend_estimators.mean_intensity_blueward[:, 0] - before_blueward,
        expected_deposition / (blend_packet.nu * abs(gradient)),
        rtol=1.0e-10,
        atol=0.0,
    )


@pytest.mark.parametrize(
    "frequency_offset", [0.0, 9.75], ids=["blend", "close-lines"]
)
def test_nonhomologous_cursor_initialization_includes_start_band(
    line_search_packet: RPacket,
    frequency_offset: float,
) -> None:
    """Claim: Cursor initialization and range admission include the close band.

    Regime: A 1e15 Hz lab frequency and 0.95e15 Hz comoving frequency, with
    exact twins or lines offset by 9.75 Hz inside the 10 Hz starting band.
    Verification: Exact cursor IDs include both starting candidates. The
    padded redward range admits the candidate 9.75 Hz above the start that
    comoving-frequency padding would exclude. This verifies bookkeeping.
    """
    frequency = 0.95e15
    frequencies = np.array(
        [
            frequency + 1.0e5,
            frequency + frequency_offset,
            frequency - frequency_offset,
            frequency - 1.0e5,
        ]
    )

    update_line_ids(line_search_packet, frequencies, frequency)

    assert line_search_packet.next_line_id == 1
    assert line_search_packet.prev_line_id == 2

    # At v_parallel = 0.05c the CMF-scaled 9.5 Hz padding is insufficient;
    # the ordinary 9.75 Hz offset belongs to the 10 Hz lab-frame close band.
    start_line_id, stop_line_id, line_id_step = get_line_id_range_nonhomologous(
        frequencies, frequency, frequency - 1.0e4, 10.0
    )
    assert (start_line_id, stop_line_id, line_id_step) == (1, 3, 1)


@pytest.mark.parametrize(
    ("active_line_id_step", "expected_ids"),
    [(1, (2, 1)), (-1, (1, 0))],
)
def test_nonhomologous_cursor_refresh_preserves_active_direction(
    line_search_packet: RPacket,
    active_line_id_step: int,
    expected_ids: tuple[int, int],
) -> None:
    """Claim: Refreshing the opposite cursor preserves the active exclusion.

    Regime: Three descending line frequencies with hand-seeded cursors that
    exclude the central line, tested for each active search direction.
    Verification: Compare the exact cursor pair with independently prescribed
    IDs: the active cursor stays fixed and the opposite one includes the
    central line. This verifies cursor bookkeeping.
    """
    frequencies = np.array([1.001e15, 1.0e15, 0.999e15])
    line_search_packet.next_line_id = 2
    line_search_packet.prev_line_id = 0

    update_line_ids(
        line_search_packet, frequencies, 1.0e15, active_line_id_step
    )

    assert (
        line_search_packet.next_line_id,
        line_search_packet.prev_line_id,
    ) == expected_ids


def test_opaque_line_interacts_in_shell_without_radial_velocity_gradient() -> None:
    """Claim: An opaque line can interact when the radial velocity gradient is zero.

    Regime: An oblique ray in a constant-velocity shell has finite projected
    gradient G = (v/r)(1 - mu**2), despite infinite radial Sobolev depth.
    Verification: Construct the opacity state with an explicit optical-depth
    coefficient, convert it through OpacityState.to_numba, and require a LINE
    event. This protects coefficient transfer against regression #3801.
    """
    velocity, radius, mu = 1.0e9, 7.5e14, 0.3
    geometry = NumbaRadial1DGeometry(
        np.array([7.0e14, 8.0e14]),
        np.array([8.0e14, 3.0e16]),
        np.array([velocity, velocity]),  # shell 0 has dv/dr = 0
        np.array([velocity, 2.0e9]),
    )
    packet_frequency = 4.0e14
    comoving_frequency = packet_frequency * (
        1.0 - velocity * mu / C_SPEED_OF_LIGHT
    )
    line_list_nu = pd.Series([comoving_frequency * (1.0 - 1.0e-7), 1.0e9])
    # Line 0 in shell 0: radial tau_sobolev = coefficient / 0 = inf.
    coefficient = pd.DataFrame([[1.0e10, 0.0], [0.0, 0.0]])
    tau_sobolev = pd.DataFrame([[np.inf, 0.0], [0.0, 0.0]])
    opacity_state = OpacityState(
        pd.Series([0.0, 0.0]),
        np.array([1.0e4, 1.0e4]),
        line_list_nu,
        tau_sobolev,
        None,
        sobolev_optical_depth_coefficient=coefficient,
    ).to_numba(None, "scatter")
    packet = RPacket(radius, mu, packet_frequency, 1.0, 42, 0)
    update_line_ids(packet, opacity_state.line_list_nu, comoving_frequency)

    _, interaction_type, _ = nonhomologous_trace_packet(
        packet,
        geometry,
        opacity_state,
        EstimatorsLine(np.zeros((2, 2)), np.zeros((2, 2))),
        0.0,  # electron opacity
        False,  # enable full relativity?
        False,  # disable line scattering?
    )

    assert interaction_type == InteractionType.LINE


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
    current_shell_id: int,
    delta_shell: int,
    no_of_shells: int,
    static_packet: RPacket,
) -> None:
    """Claim: Crossing beyond the outermost shell marks the packet as emitted.

    Regime: Prescribed shell IDs and outward steps whose destination is at
    or beyond the number of shells.
    Verification: The packet status is exactly EMITTED after the crossing.
    This verifies boundary bookkeeping.
    """
    packet = static_packet
    packet.current_shell_id = current_shell_id
    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)
    assert packet.status == PacketStatus.EMITTED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -133, 132), (132, -133, 133), (132, -1e9, 133)],
)
def test_nonhomologous_move_packet_across_shell_boundary_reabsorbed(
    current_shell_id: int,
    delta_shell: float,
    no_of_shells: int,
    static_packet: RPacket,
) -> None:
    """Claim: Crossing inside the innermost shell marks the packet as reabsorbed.

    Regime: Prescribed shell IDs and inward steps with a negative destination
    shell ID.
    Verification: The packet status is exactly REABSORBED after the crossing.
    This verifies boundary bookkeeping.
    """
    packet = static_packet
    packet.current_shell_id = current_shell_id
    move_packet_across_shell_boundary(packet, delta_shell, no_of_shells)
    assert packet.status == PacketStatus.REABSORBED


@pytest.mark.parametrize(
    ["current_shell_id", "delta_shell", "no_of_shells"],
    [(132, -1, 199), (132, 0, 132), (132, 20, 154)],
)
def test_nonhomologous_move_packet_across_shell_boundary_increment(
    current_shell_id: int,
    delta_shell: int,
    no_of_shells: int,
    static_packet: RPacket,
) -> None:
    """Claim: An internal crossing updates the shell ID by the requested step.

    Regime: Prescribed shell IDs and steps whose destination remains inside
    the shell array, including a zero step.
    Verification: The final shell ID equals the initial ID plus the step.
    This verifies boundary bookkeeping.
    """
    packet = static_packet
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
