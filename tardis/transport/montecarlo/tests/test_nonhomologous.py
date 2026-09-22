import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.opacities.macro_atom.macroatom_solver import LegacyMacroAtomSolver
from tardis.opacities.tau_sobolev import (
    calculate_sobolev_line_opacity as classic_calculate_sobolev_line_opacity,
)
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_line,
)
from tardis.transport.montecarlo.configuration.constants import (
    C_SPEED_OF_LIGHT,
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
from tardis.transport.montecarlo.modes.nonhomologous.tau_sobolev import (
    calculate_beta_sobolev,
    calculate_beta_sobolev_directional,
)
from tardis.transport.montecarlo.modes.nonhomologous.tau_sobolev import (
    calculate_sobolev_line_opacity as nonhomologous_calculate_sobolev_line_opacity,
)
from tardis.transport.montecarlo.nonhomologous_grid import (
    depressed_quartic,
    solve_resonance_quartic,
)
from tardis.transport.montecarlo.packets.movement import (
    move_packet_across_shell_boundary,
)
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
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
    assert calculate_distance_line.targetoptions["fastmath"] is False


def test_calculate_distance_line_preserves_near_line_center_root() -> None:
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

    distance = calculate_distance_line(packet, geometry, line_frequency)

    npt.assert_allclose(
        distance,
        0.14167578672128847 * shell_width,
        rtol=1.0e-9,
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

    npt.assert_array_equal(actual.to_numpy(), [[1.0], [0.0]])


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
