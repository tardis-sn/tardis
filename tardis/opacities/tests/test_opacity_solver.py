import pytest
import numpy as np
import numpy.testing as npt
import pandas.testing as pdt
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.macro_atom.macroatom_solver import (
    LegacyMacroAtomSolver,
    BoundBoundMacroAtomSolver,
)
from tardis.opacities.opacity_state import OpacityState
from tardis.opacities.tau_sobolev import calculate_sobolev_line_opacity


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
def test_opacity_solver(
    nb_simulation_verysimple, line_interaction_type, disable_line_scattering
):
    legacy_plasma = nb_simulation_verysimple.plasma

    opacity_solver = OpacitySolver(
        line_interaction_type=line_interaction_type,
        disable_line_scattering=disable_line_scattering,
    )
    actual = opacity_solver.legacy_solve(legacy_plasma)

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


@pytest.mark.parametrize(
    "line_interaction_type,disable_line_scattering",
    [
        ("macroatom", False),
        ("macroatom", True),
        ("downbranch", False),
        ("downbranch", True),
    ],
)
def test_new_macro_atom_solver(
    nb_simulation_verysimple,
    line_interaction_type,
    disable_line_scattering,
    regression_data,
):
    legacy_plasma = nb_simulation_verysimple.plasma

    opacity_solver = OpacitySolver(
        line_interaction_type=line_interaction_type,
        disable_line_scattering=disable_line_scattering,
    )
    legacy_opacity = opacity_solver.legacy_solve(legacy_plasma)

    legacy_macro_atom_state = LegacyMacroAtomSolver().solve(
        legacy_plasma.j_blues,
        legacy_plasma.atomic_data,
        legacy_opacity.tau_sobolev,
        legacy_plasma.stimulated_emission_factor,
        beta_sobolev=legacy_opacity.beta_sobolev,
    )

    macro_atom_state = BoundBoundMacroAtomSolver(
        legacy_plasma.atomic_data.levels, legacy_plasma.atomic_data.lines
    ).solve(
        legacy_plasma.j_blues,
        legacy_opacity.beta_sobolev,
        legacy_plasma.stimulated_emission_factor,
    )

    macro_atom_sorted_to_legacy = macro_atom_state.sort_to_legacy(
        legacy_macro_atom_state, legacy_plasma.atomic_data.lines
    )

    npt.assert_allclose(
        macro_atom_sorted_to_legacy.transition_probabilities.values,
        legacy_macro_atom_state.transition_probabilities.values,
        rtol=1e-7,
    )
    npt.assert_array_equal(
        macro_atom_sorted_to_legacy.transition_metadata.transition_type.values,
        legacy_plasma.atomic_data.macro_atom_data["transition_type"].values,
    )
    npt.assert_array_equal(
        macro_atom_sorted_to_legacy.transition_metadata.transition_line_idx.values,
        legacy_plasma.atomic_data.macro_atom_data["lines_idx"].values,
    )

    macro_atom_recreated = macro_atom_state.recreate_legacy_macro_atom_state(
        legacy_macro_atom_state, legacy_plasma.atomic_data.lines
    )

    npt.assert_allclose(
        macro_atom_recreated.transition_probabilities.values,
        macro_atom_sorted_to_legacy.transition_probabilities.values,
        rtol=1e-7,
    )

    npt.assert_array_equal(
        macro_atom_recreated.line2macro_level_upper,
        legacy_plasma.atomic_data.lines_upper2macro_reference_idx,
    )

    npt.assert_array_equal(
        macro_atom_recreated.macro_block_references.values,
        legacy_plasma.atomic_data.macro_atom_references[
            "block_references"
        ].values,
    )

    npt.assert_array_equal(
        macro_atom_recreated.destination_level_id.values,
        legacy_plasma.atomic_data.macro_atom_data[
            "destination_level_idx"
        ].values,
    )

    regression_new_macro_atom_data = regression_data.sync_dataframe(
        macro_atom_state.transition_probabilities
    )
    pdt.assert_frame_equal(
        regression_new_macro_atom_data,
        macro_atom_state.transition_probabilities,
    )
