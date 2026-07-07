import astropy.units as u
import numpy as np
import numpy.testing as npt
import pandas.testing as pdt
import pytest
from numpy.testing import assert_almost_equal

from tardis.opacities.macro_atom.macroatom_solver import LegacyMacroAtomSolver
from tardis.opacities.tau_sobolev import (
    calculate_sobolev_line_opacity as classic_calculate_sobolev_line_opacity,
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
from tardis.transport.montecarlo.packets.radiative_movement import (
    increment_packet_cell_index,
)
from tardis.transport.montecarlo.modes.nonhomologous.tau_sobolev import (
    calculate_beta_sobolev,
)
from tardis.transport.montecarlo.modes.nonhomologous.tau_sobolev import (
    calculate_sobolev_line_opacity as nonhomologous_calculate_sobolev_line_opacity,
)
from tardis.transport.montecarlo.nonhomologous_grid import depressed_quartic
from tardis.transport.montecarlo.packets.radiative_packet import (
    PacketStatus,
    RPacket,
)


@pytest.mark.parametrize(
    ["A", "B", "C", "D", "E", "expected_roots"],
    [
        # x^4 - 14x^3 + 71x^2 - 154x + 120 = 0
        # roots 2, 3, 4, 5
        (1.0, -14.0, 71.0, -154.0, 120.0, [5.0, 4.0, 3.0, 2.0]),
        # x^4 - x^3 = 0
        # x^3(x + 1) = 0; only real root other than 0 is 1.0
        (1.0, -1.0, 0.0, 0.0, 0.0, [1.0, 0.0, 0.0, 0.0]),
        # x^4 - 10x^3 + 35x^2 - 50x + 24 = 0
        # roots 1, 2, 3, 4
        (1.0, -10.0, 35.0, -50.0, 24.0, [4.0, 3.0, 2.0, 1.0]),
    ],
)
def test_depressed_quartic(A, B, C, D, E, expected_roots):
    """
    Standalone unit test to check accurate calculation of expected roots.
    Not a regression test.
    """
    roots = depressed_quartic(A, B, C, D, E)
    assert_almost_equal(roots, expected_roots)


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
    increment_packet_cell_index(packet, delta_shell, no_of_shells)
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
    increment_packet_cell_index(packet, delta_shell, no_of_shells)
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
    increment_packet_cell_index(packet, delta_shell, no_of_shells)
    assert packet.current_shell_id == current_shell_id + delta_shell


def test_nonhomologous_thomson_scatter(
    packet, verysimple_numba_nonhomologous_geometry
):
    """
    Analogous to
    tardis/transport/montecarlo/tests/test_interaction.py@test_thomson_scatter
    """
    init_mu = packet.mu
    init_nu = packet.nu
    init_energy = packet.energy

    nonhomologous_thomson_scatter(
        packet, verysimple_numba_nonhomologous_geometry, False
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
    verysimple_numba_nonhomologous_geometry,
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
        verysimple_numba_nonhomologous_geometry,
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
    verysimple_numba_nonhomologous_geometry,
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
        verysimple_numba_nonhomologous_geometry,
        verysimple_opacity_state,
        full_relativity,
    )

    assert packet.next_line_id == emission_line_id + 1
    npt.assert_almost_equal(packet.mu, expected["mu"])
    npt.assert_almost_equal(packet.energy, expected["energy"])


def test_nonhomologous_calculate_sobolev_line_opacity(
    nb_simulation_verysimple,
    verysimple_numba_nonhomologous_geometry,
    regression_data,
):
    """
    Analogous to
    tardis/opacities/tests/test_tau_sobolev.py@test_calculate_sobolev_line_opacity
    """
    legacy_plasma = nb_simulation_verysimple.plasma
    velocity_gradient = (
        verysimple_numba_nonhomologous_geometry.velocity_gradient
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
    verysimple_numba_nonhomologous_geometry,
    line_interaction_type,
    disable_line_scattering,
):
    """
    Analogous to
    tardis/opacities/tests/test_opacity_solver.py@test_opacity_solver
    """
    legacy_plasma = nb_simulation_verysimple.plasma
    velocity_gradient = (
        verysimple_numba_nonhomologous_geometry.velocity_gradient
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
            macro_atom_state.macro_block_references,
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
