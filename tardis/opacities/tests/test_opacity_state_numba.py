import pytest
from tardis.opacities.opacity_state import opacity_state_to_numba
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.macro_atom.macroatom_solver import MacroAtomSolver
import numpy.testing as npt
import numpy as np


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
    opacity_state = opacity_solver.solve(legacy_plasma)
    if line_interaction_type in ("downbranch", "macroatom"):
        macro_atom_state = MacroAtomSolver().solve(
            legacy_plasma,
            legacy_plasma.atomic_data,
            opacity_state.tau_sobolev,
            legacy_plasma.stimulated_emission_factor,
        )
    else:
        macro_atom_state = None
    actual = opacity_state_to_numba(
        opacity_state, macro_atom_state, line_interaction_type
    )

    if sliced:
        index = slice(2, 5)
        actual = actual[index]
    else:
        index = ...

    npt.assert_allclose(
        actual.electron_density, legacy_plasma.electron_densities.values[index]
    )
    npt.assert_allclose(
        actual.line_list_nu, legacy_plasma.atomic_data.lines.nu.values
    )
    npt.assert_allclose(
        actual.tau_sobolev, legacy_plasma.tau_sobolevs.values[:, index]
    )
    if line_interaction_type == "scatter":
        empty = np.zeros(1, dtype=np.int64)
        npt.assert_allclose(
            actual.transition_probabilities, np.zeros((1, 1), dtype=np.float64)
        )
        npt.assert_allclose(actual.line2macro_level_upper, empty)
        npt.assert_allclose(actual.macro_block_references, empty)
        npt.assert_allclose(actual.transition_type, empty)
        npt.assert_allclose(actual.destination_level_id, empty)
        npt.assert_allclose(actual.transition_line_id, empty)
    else:
        npt.assert_allclose(
            actual.transition_probabilities,
            legacy_plasma.transition_probabilities.values[:, index],
        )
        npt.assert_allclose(
            actual.line2macro_level_upper,
            legacy_plasma.atomic_data.lines_upper2macro_reference_idx,
        )
        npt.assert_allclose(
            actual.macro_block_references,
            legacy_plasma.atomic_data.macro_atom_references[
                "block_references"
            ].values,
        )
        npt.assert_allclose(
            actual.transition_type,
            legacy_plasma.atomic_data.macro_atom_data["transition_type"].values,
        )
        npt.assert_allclose(
            actual.destination_level_id,
            legacy_plasma.atomic_data.macro_atom_data[
                "destination_level_idx"
            ].values,
        )
        npt.assert_allclose(
            actual.transition_line_id,
            legacy_plasma.atomic_data.macro_atom_data["lines_idx"].values,
        )
