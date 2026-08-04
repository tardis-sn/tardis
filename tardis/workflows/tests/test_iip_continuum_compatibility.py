from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis.iip_plasma.continuum.base_continuum import BaseContinuum
from tardis.iip_plasma.continuum.base_continuum_data import ContinuumData
from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray
from tardis.io.atom_data.parse_atom_data import parse_atom_data
from tardis.io.configuration.config_reader import Configuration
from tardis.model import SimulationState
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.continuum_state import (
    EquilibriumContinuumState,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticCorrectedPhotoionizationCoeffSolver,
    CollisionalIonizationSeaton,
    SpontaneousRecombinationCoeffSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow


@dataclass(frozen=True)
class ContinuumComparisonState:
    plasma: LegacyPlasmaArray
    continuum: BaseContinuum
    photoionization_data: pd.DataFrame
    photoionization_index: pd.MultiIndex
    upper_ion_index: pd.MultiIndex
    radiation_field: DilutePlanckianRadiationField
    electron_temperature: u.Quantity
    electron_distribution: ThermalElectronEnergyDistribution
    level_to_ion_population_factor: pd.DataFrame


@dataclass(frozen=True)
class CollisionalBoundRates:
    excitation: pd.DataFrame
    deexcitation: pd.DataFrame
    excitation_index: pd.MultiIndex
    deexcitation_index: pd.MultiIndex


@pytest.fixture(scope="module")
def continuum_comparison_state(
    tardis_regression_path: Path,
) -> ContinuumComparisonState:
    """Build the restored IIP continuum used as an independent oracle."""
    comparison_config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )
    comparison_config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    comparison_config.plasma.nlte.species = [(1, 0)]

    atom_data = parse_atom_data(comparison_config)
    simulation_state = SimulationState.from_config(
        comparison_config,
        atom_data=atom_data,
    )
    atom_data.prepare_atom_data(
        simulation_state.abundance.index,
        "macroatom",
        [(1, 0)],
        [(1, 0)],
    )
    atom_data.continuum_data = ContinuumData(
        atom_data,
        selected_continuum_species=[(1, 0)],
    )
    atom_data.continuum_data.photoionization_data.loc[(1, 0, 0), "x_sect"] *= (
        0.0
    )
    atom_data.yg_data.columns = list(atom_data.collision_data_temperatures)
    atom_data.nlte_data._init_indices()
    atom_data.has_collision_data = False

    elemental_number_density = (
        simulation_state.calculate_elemental_number_density(
            atom_data.atom_data.mass
        )
    )
    plasma = LegacyPlasmaArray(
        elemental_number_density,
        atom_data,
        comparison_config.supernova.time_explosion.to("s"),
        nlte_config=comparison_config.plasma.nlte,
        delta_treatment=None,
        ionization_mode="nlte",
        excitation_mode=comparison_config.plasma.excitation,
        line_interaction_type=comparison_config.plasma.line_interaction_type,
        link_t_rad_t_electron=(
            comparison_config.plasma.link_t_rad_t_electron
            * np.ones(simulation_state.geometry.no_of_shells_active)
        ),
        helium_treatment="none",
        heating_rate_data_file=None,
        v_inner=None,
        v_outer=None,
        continuum_treatment=True,
    )
    t_radiative, dilution_factor = TypeIIPWorkflow.initialize_radiation_field(
        simulation_state.geometry,
        elemental_number_density,
        comparison_config.plasma.initial_t_inner,
        simulation_state.dilution_factor,
    )
    radiation_field = DilutePlanckianRadiationField(
        temperature=t_radiative,
        dilution_factor=dilution_factor,
    )
    j_blues = radiation_field.calculate_mean_intensity(
        plasma.atomic_data.lines.nu.values
    )
    plasma.update_radiationfield(
        t_radiative,
        dilution_factor,
        pd.DataFrame(j_blues, index=plasma.atomic_data.lines.index),
        comparison_config.plasma.nlte,
        initialize_nlte=True,
        n_e_convergence_threshold=0.05,
    )
    continuum = BaseContinuum(
        plasma_array=plasma,
        atom_data=atom_data,
        ws=dilution_factor,
        radiative_transition_probabilities=plasma.transition_probabilities,
        estimators=None,
    )

    photoionization_data = continuum.input.photoionization_data
    photoionization_index = photoionization_data.index.unique()
    upper_ion_index = pd.MultiIndex.from_arrays(
        [
            photoionization_index.get_level_values("atomic_number"),
            photoionization_index.get_level_values("ion_number") + 1,
        ],
        names=["atomic_number", "ion_number"],
    ).unique()
    electron_temperature = np.asarray(plasma.t_electrons) * u.K
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        electron_temperature,
        plasma.electron_densities.to_numpy() / u.cm**3,
    )
    level_to_ion_population_factor = plasma.lte_level_number_density.loc[
        photoionization_index
    ].divide(
        plasma.lte_ion_number_density.loc[upper_ion_index].to_numpy()
        * plasma.electron_densities.to_numpy(),
    )

    return ContinuumComparisonState(
        plasma=plasma,
        continuum=continuum,
        photoionization_data=photoionization_data,
        photoionization_index=photoionization_index,
        upper_ion_index=upper_ion_index,
        radiation_field=radiation_field,
        electron_temperature=electron_temperature,
        electron_distribution=electron_distribution,
        level_to_ion_population_factor=level_to_ion_population_factor,
    )


@pytest.fixture(scope="module")
def assembled_equilibrium_continuum_state(
    continuum_comparison_state: ContinuumComparisonState,
) -> EquilibriumContinuumState:
    """Assemble the replacement continuum state from the oracle plasma."""
    state = continuum_comparison_state
    legacy_atomic_data = state.plasma.atomic_data
    atomic_data = SimpleNamespace(
        photoionization_data=state.photoionization_data,
        levels=legacy_atomic_data.levels,
        lines=legacy_atomic_data.lines,
        collision_data_temperatures=(
            legacy_atomic_data.collision_data_temperatures
        ),
        yg_data=legacy_atomic_data.yg_data,
        ionization_data=legacy_atomic_data.ionization_data,
    )
    plasma = SimpleNamespace(
        atomic_data=atomic_data,
        t_electrons=state.plasma.t_electrons,
        electron_densities=state.plasma.electron_densities,
        dilute_planckian_radiation_field=state.radiation_field,
        lte_level_number_density=state.plasma.lte_level_number_density,
        level_number_density=state.plasma.level_number_density,
        lte_ion_number_density=state.plasma.lte_ion_number_density,
        ion_number_density=state.plasma.ion_number_density,
    )
    return EquilibriumContinuumState.from_plasma(plasma)


@pytest.fixture(scope="module")
def collisional_bound_rates(
    continuum_comparison_state: ContinuumComparisonState,
) -> CollisionalBoundRates:
    state = continuum_comparison_state
    plasma = state.plasma
    collisional_rates = ThermalCollisionalRateSolver(
        plasma.atomic_data.levels,
        plasma.atomic_data.lines,
        plasma.atomic_data.collision_data_temperatures,
        plasma.atomic_data.yg_data,
        collision_strengths_type="cmfgen",
    ).solve(state.electron_temperature)
    source_levels = collisional_rates.index.get_level_values(
        "level_number_source"
    )
    destination_levels = collisional_rates.index.get_level_values(
        "level_number_destination"
    )
    excitation = collisional_rates[source_levels < destination_levels]
    deexcitation = collisional_rates[source_levels > destination_levels]
    excitation_index = excitation.index.droplevel(
        ["ion_number_source", "ion_number_destination"]
    ).rename(
        {
            "level_number_source": "level_number_lower",
            "level_number_destination": "level_number_upper",
        }
    )
    deexcitation_index = (
        deexcitation.index.swaplevel(
            "level_number_source", "level_number_destination"
        )
        .droplevel(["ion_number_source", "ion_number_destination"])
        .rename(
            {
                "level_number_destination": "level_number_lower",
                "level_number_source": "level_number_upper",
            }
        )
    )
    return CollisionalBoundRates(
        excitation,
        deexcitation,
        excitation_index,
        deexcitation_index,
    )


@pytest.fixture(scope="module")
def collisional_ionization_rate(
    continuum_comparison_state: ContinuumComparisonState,
) -> pd.DataFrame:
    return CollisionalIonizationSeaton(
        continuum_comparison_state.photoionization_data
    ).solve(continuum_comparison_state.electron_temperature)


@pytest.fixture(scope="module")
def equilibrium_cooling_channels(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_bound_rates: CollisionalBoundRates,
    collisional_ionization_rate: pd.DataFrame,
) -> npt.NDArray[np.float64]:
    state = continuum_comparison_state
    plasma = state.plasma
    bound_free_cooling = BoundFreeThermalRates(
        state.photoionization_data
    ).solve(
        plasma.level_number_density,
        plasma.ion_number_density,
        state.electron_distribution,
        state.level_to_ion_population_factor,
        state.radiation_field,
    )[1]
    free_free_cooling = FreeFreeThermalRates().solve(
        pd.Series(0.0, index=plasma.electron_densities.index),
        state.electron_distribution,
        plasma.ion_number_density,
    )[1]
    collisional_ionization_cooling = CollisionalIonizationThermalRates(
        state.photoionization_data
    ).solve(
        state.electron_distribution.number_density,
        plasma.ion_number_density,
        plasma.level_number_density,
        collisional_ionization_rate,
        state.level_to_ion_population_factor,
    )[1]
    collisional_bound_cooling = CollisionalBoundThermalRates(
        plasma.atomic_data.lines.loc[collisional_bound_rates.excitation_index]
    ).solve(
        state.electron_distribution.number_density,
        collisional_bound_rates.deexcitation.set_axis(
            collisional_bound_rates.deexcitation_index
        ),
        collisional_bound_rates.excitation.set_axis(
            collisional_bound_rates.excitation_index
        ),
        plasma.level_number_density,
    )[1]

    return np.vstack(
        [
            collisional_bound_cooling.to_numpy(),
            collisional_ionization_cooling.to_numpy(),
            bound_free_cooling.to_numpy(),
            free_free_cooling.to_numpy(),
        ]
    )


def test_radiative_ionization_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
) -> None:
    state = continuum_comparison_state
    radiative_ionization_rate = AnalyticCorrectedPhotoionizationCoeffSolver(
        state.photoionization_data
    ).solve(
        state.radiation_field,
        state.electron_temperature,
        state.plasma.lte_level_number_density.loc[state.photoionization_index],
        state.plasma.level_number_density.loc[state.photoionization_index],
        state.plasma.lte_ion_number_density.loc[state.upper_ion_index],
        state.plasma.ion_number_density.loc[state.upper_ion_index],
    )
    pd.testing.assert_index_equal(
        radiative_ionization_rate.index,
        state.continuum.radiative_ionization.rate_coefficient.index,
    )
    np.testing.assert_allclose(
        radiative_ionization_rate.to_numpy(),
        state.continuum.radiative_ionization.rate_coefficient.to_numpy(),
        rtol=2e-6,
        atol=0.0,
    )


def test_radiative_recombination_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
) -> None:
    state = continuum_comparison_state
    radiative_recombination_rate = (
        SpontaneousRecombinationCoeffSolver(state.photoionization_data).solve(
            state.electron_temperature
        )
        * state.level_to_ion_population_factor
    )
    pd.testing.assert_index_equal(
        radiative_recombination_rate.index,
        state.continuum.radiative_recombination.rate_coefficient.index,
    )
    np.testing.assert_allclose(
        radiative_recombination_rate.to_numpy(),
        state.continuum.radiative_recombination.rate_coefficient.to_numpy(),
        rtol=2e-4,
        atol=0.0,
    )


def test_collisional_excitation_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_bound_rates: CollisionalBoundRates,
) -> None:
    np.testing.assert_allclose(
        collisional_bound_rates.excitation.to_numpy(),
        continuum_comparison_state.plasma.coll_exc_coeff.loc[
            collisional_bound_rates.excitation_index
        ].to_numpy(),
        rtol=2e-5,
        atol=0.0,
    )


def test_collisional_deexcitation_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_bound_rates: CollisionalBoundRates,
) -> None:
    np.testing.assert_allclose(
        collisional_bound_rates.deexcitation.to_numpy(),
        continuum_comparison_state.plasma.coll_deexc_coeff.loc[
            collisional_bound_rates.deexcitation_index
        ].to_numpy(),
        rtol=2e-5,
        atol=0.0,
    )


def test_collisional_ionization_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_ionization_rate: pd.DataFrame,
) -> None:
    pd.testing.assert_frame_equal(
        collisional_ionization_rate,
        continuum_comparison_state.plasma.coll_ion_coeff.loc[
            collisional_ionization_rate.index
        ],
        check_names=False,
        check_column_type=False,
    )


def test_collisional_recombination_rates_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_ionization_rate: pd.DataFrame,
) -> None:
    collisional_recombination_rate = (
        collisional_ionization_rate
        * continuum_comparison_state.level_to_ion_population_factor
    )
    np.testing.assert_allclose(
        collisional_recombination_rate.to_numpy(),
        continuum_comparison_state.plasma.coll_recomb_coeff.loc[
            collisional_recombination_rate.index
        ].to_numpy(),
        rtol=2e-5,
        atol=0.0,
    )


def test_assembled_radiative_recombination_rate_matches_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    assembled_equilibrium_continuum_state: EquilibriumContinuumState,
) -> None:
    np.testing.assert_allclose(
        assembled_equilibrium_continuum_state.radiative_recombination_rate,
        continuum_comparison_state.continuum.radiative_recombination.rate_coefficient,
        rtol=2e-4,
        atol=0.0,
    )


def test_assembled_collisional_recombination_rate_matches_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    assembled_equilibrium_continuum_state: EquilibriumContinuumState,
) -> None:
    actual = (
        assembled_equilibrium_continuum_state.collisional_recombination_rate
    )
    expected = continuum_comparison_state.plasma.coll_recomb_coeff
    pd.testing.assert_index_equal(actual.index, expected.index)
    np.testing.assert_allclose(
        actual,
        expected,
        rtol=2e-5,
        atol=0.0,
    )


def test_assembled_collisional_deexcitation_rate_matches_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_bound_rates: CollisionalBoundRates,
    assembled_equilibrium_continuum_state: EquilibriumContinuumState,
) -> None:
    actual = assembled_equilibrium_continuum_state.collisional_deexcitation_rate
    expected_index = collisional_bound_rates.deexcitation_index
    expected = continuum_comparison_state.plasma.coll_deexc_coeff.loc[
        expected_index
    ]
    pd.testing.assert_index_equal(actual.index, expected_index)
    np.testing.assert_allclose(actual, expected, rtol=2e-5, atol=0.0)


def test_assembled_collision_energy_gaps_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    collisional_bound_rates: CollisionalBoundRates,
    assembled_equilibrium_continuum_state: EquilibriumContinuumState,
) -> None:
    actual = assembled_equilibrium_continuum_state.delta_E_yg
    expected = continuum_comparison_state.plasma.delta_E_yg.loc[
        collisional_bound_rates.deexcitation_index
    ]
    pd.testing.assert_series_equal(actual, expected, check_names=False)
    assert np.all(actual > 0.0)


def test_assembled_free_bound_cooling_probability_matches_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    assembled_equilibrium_continuum_state: EquilibriumContinuumState,
) -> None:
    expected = continuum_comparison_state.continuum.cooling_rates
    np.testing.assert_allclose(
        assembled_equilibrium_continuum_state.radiative_recombination_cooling_probability,
        expected.radiative_recombination_probability,
        rtol=3e-4,
        atol=0.0,
    )


def test_assembled_free_bound_cooling_distribution_matches_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    assembled_equilibrium_continuum_state: EquilibriumContinuumState,
) -> None:
    expected = continuum_comparison_state.continuum.cooling_rates
    np.testing.assert_allclose(
        assembled_equilibrium_continuum_state.radiative_recombination_cooling_array,
        expected.radiative_recombination.probabilities_array,
        rtol=3e-4,
        atol=0.0,
    )


def test_cooling_channel_totals_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    equilibrium_cooling_channels: npt.NDArray[np.float64],
) -> None:
    cooling_rates = continuum_comparison_state.continuum.cooling_rates
    iip_cooling_channels = np.vstack(
        [
            cooling_rates.collisional_excitation_total,
            cooling_rates.collisional_ionization_total,
            cooling_rates.radiative_recombination_total,
            cooling_rates.free_free_total,
        ]
    )
    np.testing.assert_allclose(
        equilibrium_cooling_channels,
        iip_cooling_channels,
        rtol=3e-4,
        atol=0.0,
    )


def test_cooling_channel_probabilities_match_iip_continuum(
    continuum_comparison_state: ContinuumComparisonState,
    equilibrium_cooling_channels: npt.NDArray[np.float64],
) -> None:
    cooling_rates = continuum_comparison_state.continuum.cooling_rates
    np.testing.assert_allclose(
        equilibrium_cooling_channels / equilibrium_cooling_channels.sum(axis=0),
        np.vstack(
            [
                cooling_rates.collisional_excitation_probability,
                cooling_rates.collisional_ionization_probability,
                cooling_rates.radiative_recombination_probability,
                cooling_rates.free_free_probability,
            ]
        ),
        rtol=3e-4,
        atol=0.0,
    )


@pytest.mark.parametrize(
    "process_name",
    [
        "collisional_excitation",
        "collisional_ionization",
        "radiative_recombination",
    ],
)
def test_iip_process_probabilities_normalize_per_shell(
    continuum_comparison_state: ContinuumComparisonState,
    process_name: str,
) -> None:
    cooling_channel = getattr(
        continuum_comparison_state.continuum.cooling_rates,
        process_name,
    )
    assert cooling_channel.probabilities_array.shape == (
        len(continuum_comparison_state.plasma.t_electrons),
        len(cooling_channel.references),
    )
    np.testing.assert_allclose(
        cooling_channel.probabilities_array.sum(axis=1),
        1.0,
        rtol=1e-12,
        atol=0.0,
    )
