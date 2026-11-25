import numpy as np
import numpy.testing as npt
import pytest

from tardis.opacities.continuum.continuum_state_numba import (
    ContinuumOpacityStateNumba,
)
from tardis.opacities.macro_atom.macroatom_solver import (
    BoundBoundMacroAtomSolver,
)
from tardis.opacities.opacities import chi_continuum_calculator
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.opacity_state_numba import OpacityStateNumba


@pytest.mark.parametrize(
    "line_interaction_type,sliced",
    [
        ("scatter", False),
        ("macroatom", False),
        ("macroatom", True),
        ("downbranch", False),
        ("downbranch", True),
    ],
)
def test_opacity_state_to_numba(
    nb_simulation_verysimple, line_interaction_type, sliced
):
    legacy_plasma = nb_simulation_verysimple.plasma

    opacity_solver = OpacitySolver(
        line_interaction_type=line_interaction_type,
        disable_line_scattering=False,
    )
    opacity_state = opacity_solver.legacy_solve(legacy_plasma)
    if line_interaction_type in ("downbranch", "macroatom"):
        macro_atom_state = BoundBoundMacroAtomSolver(
            legacy_plasma.atomic_data.levels, legacy_plasma.atomic_data.lines
        ).solve(
            legacy_plasma.j_blues,
            opacity_state.beta_sobolev,
            legacy_plasma.stimulated_emission_factor,
        )
    else:
        macro_atom_state = None
    actual = opacity_state.to_numba(macro_atom_state, line_interaction_type)

    if sliced:
        index = slice(2, 5)
        actual = actual[index]
    else:
        index = ...

    npt.assert_allclose(
        actual.electron_density,
        legacy_plasma.electron_densities.to_numpy()[index],
    )
    npt.assert_allclose(
        actual.line_list_nu, legacy_plasma.atomic_data.lines.nu.to_numpy()
    )
    npt.assert_allclose(
        actual.tau_sobolev, legacy_plasma.tau_sobolevs.to_numpy()[:, index]
    )
    if line_interaction_type == "scatter":
        empty = np.zeros(1, dtype=np.int64)
        npt.assert_allclose(
            actual.transition_probabilities, np.zeros((1, 1), dtype=np.float64)
        )
        npt.assert_allclose(actual.line2macro_level_upper, empty)
        npt.assert_allclose(actual.macro_block_edge_index, empty)
        npt.assert_allclose(actual.transition_type, empty)
        npt.assert_allclose(actual.destination_level_id, empty)
        npt.assert_allclose(actual.transition_line_id, empty)


def test_sliced_opacity_state_keeps_each_shell_continuum_opacity() -> None:
    """Slicing selects shells without changing any shell's continuum opacity.

    Shell k of the slice [1:3] must give the bound-free and free-free opacity
    of shell 1 + k of the full state. Every shell has its own electron
    temperature, free-free factor and bound-free table, so pairing data from
    different shells changes the opacity. Both sides run the same arithmetic
    on the same numbers, so equality is exact.
    """
    n_shells = 4
    continuum = ContinuumOpacityStateNumba(
        np.array([3.0e15]),
        np.ones((1, n_shells)),
        np.array([3.0e15]),
        np.array([4.0e15]),
        np.array([0, 2]),
        # Bound-free opacity [cm^-1] at the two tabulated frequencies, per shell.
        np.array([[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]]) * 1.0e-15,
        np.array([1.0e-18, 2.0e-18]),
        np.array([3.0e15, 4.0e15]),
        np.array([1.0, 2.0, 3.0, 4.0]),
        np.ones((1, n_shells)),
        np.zeros(1, dtype=np.int64),
        0,
        np.ones((n_shells, 1, 1)),
        np.arange(1, dtype=np.int64),
        np.arange(1, dtype=np.int64),
    )
    no_macro_atom = np.zeros(1, dtype=np.int64)
    opacity_state = OpacityStateNumba(
        np.full(n_shells, 1.0e8),
        np.array([5.0e3, 1.0e4, 2.0e4, 4.0e4]),
        np.array([1.0e15]),
        np.zeros((1, n_shells)),
        np.zeros((1, n_shells)),
        no_macro_atom,
        no_macro_atom,
        no_macro_atom,
        no_macro_atom,
        no_macro_atom,
        continuum,
    )
    sliced = opacity_state[1:3]
    nu = 3.5e15

    for shell in range(2):
        expected = chi_continuum_calculator(
            opacity_state.continuum, opacity_state, nu, 1 + shell
        )
        actual = chi_continuum_calculator(sliced.continuum, sliced, nu, shell)
        # (bound-free total, free-free)
        assert (actual[0], actual[4]) == (expected[0], expected[4])
