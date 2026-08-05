from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tardis import constants as const
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma import BasePlasma
from tardis.plasma.equilibrium.continuum_state import (
    EquilibriumContinuumState,
)
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow


@pytest.fixture(scope="module")
def continuum_workflow(tardis_regression_path: Path) -> TypeIIPWorkflow:
    """Build the real standard IIP workflow used by continuum tests."""
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )
    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    config.plasma.nlte.species = [(1, 0)]
    return TypeIIPWorkflow(config)


@pytest.fixture(scope="module")
def helium_continuum_workflow(tardis_regression_path: Path) -> TypeIIPWorkflow:
    """Build a real standard workflow with helium continuum interactions."""
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )
    config.atom_data = (
        tardis_regression_path / "atom_data" / "TestNLTE_He_Ti.h5"
    )
    config.model.abundances.pop("H", None)
    config.model.abundances["He"] = 1.0
    config.plasma.continuum_interaction.species = ["He I"]
    config.plasma.nlte.species = [(2, 0)]
    return TypeIIPWorkflow(config)


@pytest.fixture(scope="module")
def hydrogen_helium_continuum_workflow(
    tardis_regression_path: Path,
) -> TypeIIPWorkflow:
    """Build a real standard workflow with H I and He I continua."""
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )
    config.atom_data = (
        tardis_regression_path / "atom_data" / "TestNLTE_He_Ti.h5"
    )
    config.model.abundances["H"] = 0.5
    config.model.abundances["He"] = 0.5
    config.plasma.continuum_interaction.species = ["H I", "He I"]
    config.plasma.nlte.species = [(1, 0), (2, 0)]
    return TypeIIPWorkflow(config)


def test_standard_continuum_state_is_complete(
    continuum_workflow: TypeIIPWorkflow,
) -> None:
    state = continuum_workflow.continuum_state
    assert isinstance(state, EquilibriumContinuumState)
    for value in (
        state.radiative_ionization_rate,
        state.radiative_recombination_rate,
        state.collisional_excitation_rate,
        state.collisional_deexcitation_rate,
        state.collisional_ionization_rate,
        state.collisional_recombination_rate,
    ):
        assert np.all(np.isfinite(value.to_numpy()))


def test_continuum_rate_indices_are_consistent(
    continuum_workflow: TypeIIPWorkflow,
) -> None:
    state = continuum_workflow.continuum_state
    assert state.radiative_recombination_rate.index.equals(
        state.radiative_ionization_rate.index
    )
    assert state.collisional_recombination_rate.index.equals(
        state.collisional_ionization_rate.index
    )
    assert state.delta_E_yg.index.equals(
        state.collisional_deexcitation_rate.index
    )
    assert np.all(state.delta_E_yg.to_numpy() > 0.0)


def test_continuum_cooling_probabilities_are_normalized(
    continuum_workflow: TypeIIPWorkflow,
) -> None:
    """Check structured cooling contracts replacing legacy cooling fields."""
    state = continuum_workflow.continuum_state
    for probability, cooling_array in (
        (
            state.collisional_excitation_cooling_probability,
            state.collisional_excitation_cooling_array,
        ),
        (
            state.collisional_ionization_cooling_probability,
            state.collisional_ionization_cooling_array,
        ),
        (
            state.radiative_recombination_cooling_probability,
            state.radiative_recombination_cooling_array,
        ),
    ):
        assert np.all(np.isfinite(np.asarray(probability)))
        np.testing.assert_allclose(cooling_array.sum(axis=1), 1.0)
    np.testing.assert_allclose(
        state.collisional_excitation_cooling_probability
        + state.collisional_ionization_cooling_probability
        + state.radiative_recombination_cooling_probability
        + state.free_free_cooling_probability,
        1.0,
    )


def test_structured_level_populations_replace_legacy_departure_coefficient(
    continuum_workflow: TypeIIPWorkflow,
) -> None:
    """Protect population physics without exposing the legacy ``b`` field."""
    plasma = continuum_workflow.plasma_solver
    level_ion_index = plasma.level_number_density.index.droplevel(
        "level_number"
    )
    expected = pd.DataFrame(
        plasma.hydrogen_continuum_level_boltzmann_factor.to_numpy()
        / plasma.hydrogen_continuum_partition_function.loc[
            level_ion_index
        ].to_numpy()
        * plasma.ion_number_density.loc[level_ion_index].to_numpy(),
        index=plasma.level_number_density.index,
        columns=plasma.level_number_density.columns,
    )
    pd.testing.assert_frame_equal(
        plasma.level_number_density,
        expected,
        check_names=False,
        check_dtype=False,
        rtol=2e-12,
        atol=0.0,
    )


def test_continuum_opacity_state_obeys_bound_free_identity(
    continuum_workflow: TypeIIPWorkflow,
) -> None:
    workflow = continuum_workflow
    opacity = workflow.continuum_opacity_state
    photo_data = workflow.plasma_solver.photo_ion_cross_sections
    photo_levels = workflow.plasma_solver.level_number_density.loc[
        photo_data.index
    ]
    lte_photo_levels = workflow.plasma_solver.lte_level_number_density.loc[
        photo_data.index
    ]
    upper_ions = pd.MultiIndex.from_arrays(
        [
            photo_data.index.get_level_values("atomic_number"),
            photo_data.index.get_level_values("ion_number") + 1,
        ],
        names=["atomic_number", "ion_number"],
    )
    ion_ratio = workflow.plasma_solver.ion_number_density.loc[
        upper_ions
    ].divide(
        workflow.plasma_solver.lte_ion_number_density.loc[upper_ions]
    )
    boltzmann_factor = np.exp(
        -photo_data.nu.to_numpy()[:, np.newaxis]
        / workflow.plasma_solver.t_electrons[np.newaxis, :]
            * (const.h.cgs.value / const.k_B.cgs.value)
    )
    expected = (
        photo_levels
        * (
            1.0
            - lte_photo_levels.divide(photo_levels)
            .multiply(ion_ratio.to_numpy())
            .multiply(boltzmann_factor)
        )
    ).multiply(photo_data.x_sect.to_numpy(), axis=0)
    pd.testing.assert_frame_equal(
        opacity._chi_bf,
        expected,
        check_names=False,
        check_dtype=False,
        rtol=2e-12,
        atol=0.0,
    )


def test_continuum_macro_atom_state_is_structured(
    continuum_workflow: TypeIIPWorkflow,
) -> None:
    workflow = continuum_workflow
    macro_atom_state = workflow.solve_macro_atom_state()
    assert macro_atom_state.transition_probabilities.shape[0] > 0
    assert macro_atom_state.absorbing_probability_matrix is not None
    assert not hasattr(workflow.plasma_solver, "p_fb_deactivation")
    assert not hasattr(workflow.plasma_solver, "chi_bf")
    assert not hasattr(workflow.plasma_solver, "transition_probabilities")
    assert np.all(
        np.isfinite(macro_atom_state.transition_probabilities.to_numpy())
    )


def test_type_iip_workflow_uses_the_standard_plasma(
    continuum_workflow: TypeIIPWorkflow,
) -> None:
    """Keep the legacy plasma implementation out of the production workflow."""
    assert type(continuum_workflow.plasma_solver) is BasePlasma


def test_standard_workflow_supports_non_hydrogen_continuum(
    helium_continuum_workflow: TypeIIPWorkflow,
) -> None:
    """Verify a real workflow can assemble and expose helium continuum data."""
    workflow = helium_continuum_workflow
    assert set(
        workflow.plasma_solver.ion_number_density.index.get_level_values(
            "atomic_number"
        )
    ) == {2}
    assert np.all(
        np.isfinite(workflow.continuum_state.radiative_ionization_rate)
    )
    assert np.all(np.isfinite(workflow.continuum_opacity_state.chi_bf))
    macro_atom_state = workflow.solve_macro_atom_state()
    assert np.all(
        np.isfinite(macro_atom_state.transition_probabilities.to_numpy())
    )


def test_standard_workflow_supports_hydrogen_and_helium_continuum(
    hydrogen_helium_continuum_workflow: TypeIIPWorkflow,
) -> None:
    """Verify finite, conserved populations for multiple continua."""
    workflow = hydrogen_helium_continuum_workflow
    plasma = workflow.plasma_solver
    assert set(
        plasma.ion_number_density.index.get_level_values("atomic_number")
    ) == {1, 2}
    level_totals = plasma.level_number_density.groupby(
        level=["atomic_number", "ion_number"]
    ).sum()
    np.testing.assert_allclose(
        level_totals.to_numpy(),
        plasma.ion_number_density.to_numpy(),
        rtol=1e-12,
        atol=0.0,
    )
    assert np.all(np.isfinite(plasma.electron_densities.to_numpy()))
    assert np.all(np.isfinite(workflow.continuum_opacity_state.chi_bf))
    macro_atom_state = workflow.solve_macro_atom_state()
    assert np.all(
        np.isfinite(macro_atom_state.transition_probabilities.to_numpy())
    )
    free_bound_destinations = macro_atom_state.transition_metadata.loc[
        macro_atom_state.transition_metadata.transition_type
        == MacroAtomTransitionType.FB_COOLING,
        "destination",
    ]
    assert {destination[0] for destination in free_bound_destinations} == {
        1,
        2,
    }
