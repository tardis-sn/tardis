from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis.iip_plasma.continuum.base_continuum import BaseContinuum
from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.level_populations import LevelPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import RateMatrix
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
from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
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


def _max_rel_diff(actual, expected):
    actual_vals = actual.values
    expected_vals = expected.values

    relative_difference = np.abs(actual_vals - expected_vals) / np.abs(
        expected_vals
    )

    return float(np.nanmax(relative_difference))


@pytest.fixture
def iip_regression_path(tardis_regression_path):
    return tardis_regression_path / "tardis" / "workflows" / "tests"


@pytest.fixture
def ctardis_compare_config(
    tardis_regression_path: Path,
) -> Configuration:
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )

    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )

    config.plasma.nlte.species = [
        (1, 0)
    ]  # Hack to force config necessary for ctardis plasma
    return config


@pytest.fixture
def type_iip_workflow(ctardis_compare_config):
    workflow = TypeIIPWorkflow(ctardis_compare_config)
    return workflow


# identical to ctardis values
@pytest.fixture
def elemental_number_density(iip_regression_path):
    elemental_number_density = pd.read_hdf(
        iip_regression_path / "ctardis_elemental_number_density.h5",
        key="data",
    )
    elemental_number_density.columns = elemental_number_density.columns.astype(
        int
    )
    return elemental_number_density


# initial plasma setup matching ctardis
@pytest.fixture
def iip_plasma(iip_atom_data, elemental_number_density, ctardis_compare_config):
    plasma = LegacyPlasmaArray(
        elemental_number_density,
        iip_atom_data,
        ctardis_compare_config.supernova.time_explosion.to("s").value,
        nlte_config=ctardis_compare_config.plasma.nlte,
        delta_treatment=None,
        ionization_mode="nlte",
        excitation_mode="dilute-lte",
        line_interaction_type=ctardis_compare_config.plasma.line_interaction_type,
        link_t_rad_t_electron=1.0 * np.ones(24),
        # link_t_rad_t_electron=self.ws**0.25,
        helium_treatment="none",
        heating_rate_data_file=None,
        v_inner=None,
        v_outer=None,
        continuum_treatment=True,
    )

    return plasma


# "NLTE init" is the first call to update_radiationfield to set up the plasma
@pytest.fixture
def iip_plasma_nlte_init(
    iip_regression_path, iip_plasma, ctardis_compare_config
):
    j_blues_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_j_blues_ctardis_init_nlte.h5",
        key="data",
    )

    # ctardis starts with a constant rad temperature in all cells
    radiation_temp = 9984.96131287 * np.ones(24)
    dilution_factor = np.array(
        [
            0.18635244,
            0.15938095,
            0.11736085,
            0.34665656,
            0.32265696,
            0.30224056,
            0.28436446,
            0.26841929,
            0.2540108,
            0.24086562,
            0.22878441,
            0.21761613,
            0.20724285,
            0.1975702,
            0.18852112,
            0.18003167,
            0.17204798,
            0.16452412,
            0.15742053,
            0.15070279,
            0.14434073,
            0.13830767,
            0.13257993,
            0.12856901,
        ]
    )

    iip_plasma.update_radiationfield(
        radiation_temp,
        dilution_factor,
        j_blues_ctardis,
        ctardis_compare_config.plasma.nlte,
        initialize_nlte=True,
        n_e_convergence_threshold=0.05,
    )
    return iip_plasma


@pytest.fixture
def iip_plasma_after_mc(
    iip_regression_path, iip_plasma_nlte_init, ctardis_compare_config
):
    j_blues_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_j_blues_ctardis_after_mc.h5",
        key="data",
    )

    radiation_temp = np.array(
        [
            9992.27229695,
            9992.59224105,
            9983.78000964,
            9980.58614386,
            9979.83477025,
            9968.05132981,
            9957.88724805,
            9949.36369847,
            9946.8743961,
            9937.71425418,
            9934.85610192,
            9928.23880919,
            9926.40535242,
            9916.93223133,
            9912.22246589,
            9911.051763,
            9910.26097021,
            9901.72775668,
            9895.9432972,
            9891.58754489,
            9886.70685954,
            9880.93185734,
            9876.00858684,
            9872.59842944,
        ]
    )

    dilution_factor = np.array(
        [
            0.3571996,
            0.31756545,
            0.27019532,
            0.36604569,
            0.33787167,
            0.31579601,
            0.29590609,
            0.27936991,
            0.2634541,
            0.24940025,
            0.23579985,
            0.22373621,
            0.21241799,
            0.20254584,
            0.19309261,
            0.18394483,
            0.1755579,
            0.16798016,
            0.16076174,
            0.15381029,
            0.14730572,
            0.14119434,
            0.13532174,
            0.13124624,
        ]
    )

    photo_ion_estimator = pd.read_hdf(
        iip_regression_path / "ctardis_photo_ion_estimator_after_mc.h5",
        key="data",
    )

    stim_recomb_estimator = pd.read_hdf(
        iip_regression_path / "ctardis_stim_recomb_estimator_after_mc.h5",
        key="data",
    )

    bf_heating_estimator = pd.read_hdf(
        iip_regression_path / "ctardis_bf_heating_estimator_after_mc.h5",
        key="data",
    )

    stim_recomb_cooling_estimator = pd.read_hdf(
        iip_regression_path
        / "ctardis_stim_recomb_cooling_estimator_after_mc.h5",
        key="data",
    )

    ff_heating_estimator = [
        4.89135279e-24,
        4.37696370e-24,
        3.75869301e-24,
        4.97847160e-24,
        4.52158002e-24,
        4.21024499e-24,
        3.94991540e-24,
        3.72915649e-24,
        3.58902110e-24,
        3.40170224e-24,
        3.20848519e-24,
        3.03540032e-24,
        2.87314722e-24,
        2.74328938e-24,
        2.61063140e-24,
        2.50640248e-24,
        2.38164559e-24,
        2.26967531e-24,
        2.24509826e-24,
        2.12378192e-24,
        2.02063266e-24,
        1.92509873e-24,
        1.83070678e-24,
        1.77346374e-24,
    ]

    continuum_estimators = {}

    continuum_estimators["photo_ion_estimator"] = photo_ion_estimator
    continuum_estimators["stim_recomb_estimator"] = stim_recomb_estimator
    continuum_estimators["bf_heating_estimator"] = bf_heating_estimator
    continuum_estimators["stim_recomb_cooling_estimator"] = (
        stim_recomb_cooling_estimator
    )
    continuum_estimators["ff_heating_estimator"] = ff_heating_estimator

    iip_plasma_nlte_init.update_radiationfield(
        radiation_temp,
        dilution_factor,
        j_blues_ctardis,
        ctardis_compare_config.plasma.nlte,
        initialize_nlte=False,
        n_e_convergence_threshold=0.05,
        **continuum_estimators,
    )

    return iip_plasma_nlte_init


def test_iip_plasma_initialization(iip_plasma_nlte_init, iip_regression_path):
    tau_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_tau_sobolevs_init_nlte.h5",
        key="data",
    )
    beta_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_beta_sobolevs_init_nlte.h5",
        key="data",
    )
    ion_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_ion_density_init_nlte.h5",
        key="data",
    )
    level_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_level_number_density_init_nlte.h5",
        key="data",
    )
    transition_probabilities_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_transition_probabilities_init_nlte.h5",
        key="data",
    )
    electron_densities_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_electron_densities_init_nlte.h5",
        key="data",
    )
    p_fb_deactivation_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_p_fb_deactivation_init_nlte.h5",
        key="data",
    )
    chi_bf_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_chi_bf_init_nlte.h5",
        key="data",
    )
    t_electrons_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_t_electrons_init_nlte.h5", key="data"
    )

    print(
        "init transition_probabilities max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.transition_probabilities,
                transition_probabilities_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_nlte_init.transition_probabilities,
        transition_probabilities_ctardis,
        rtol=3e-8,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "init ion_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.ion_number_density,
                ion_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_nlte_init.ion_number_density,
        ion_number_density_ctardis,
        rtol=4e-8,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "init tau_sobolevs max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.tau_sobolevs,
                tau_sobolevs_ctardis,
            )
        )
    )
    # Sobolev values are stored differently between codes, so comparing raw data instead
    np.testing.assert_allclose(
        iip_plasma_nlte_init.tau_sobolevs.values,
        tau_sobolevs_ctardis.values,
        rtol=4e-8,
        atol=0,
    )

    print(
        "init beta_sobolev max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.beta_sobolev,
                beta_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_nlte_init.beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=3e-8,
        atol=0,
    )

    print(
        "init level_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.level_number_density,
                level_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_nlte_init.level_number_density,
        level_number_density_ctardis,
        rtol=4e-8,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "init electron_densities max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.electron_densities,
                electron_densities_ctardis,
            )
        )
    )
    pd.testing.assert_series_equal(
        iip_plasma_nlte_init.electron_densities,
        electron_densities_ctardis,
        rtol=2e-12,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "init t_electrons max rel diff: {:.3e}".format(
            _max_rel_diff(
                pd.DataFrame(iip_plasma_nlte_init.t_electrons),
                t_electrons_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_nlte_init.t_electrons,
        t_electrons_ctardis.values,
        rtol=2e-13,
        atol=0,
    )

    print(
        "init p_fb_deactivation max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.p_fb_deactivation,
                p_fb_deactivation_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_nlte_init.p_fb_deactivation.values,
        p_fb_deactivation_ctardis.values,
        rtol=2e-13,
        atol=0,
    )

    print(
        "init chi_bf max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_nlte_init.chi_bf,
                chi_bf_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_nlte_init.chi_bf.values,
        chi_bf_ctardis.values,
        rtol=4e-8,
        atol=0,
    )


@pytest.fixture(scope="module")
def continuum_comparison_state(
    tardis_regression_path: Path,
) -> ContinuumComparisonState:
    comparison_config = Configuration.from_yaml(
        "tardis/workflows/tests/data/ctardis_compare.yml"
    )
    comparison_config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    comparison_config.plasma.nlte.species = [(1, 0)]
    workflow = TypeIIPWorkflow(comparison_config)
    plasma = workflow.plasma_solver
    continuum = workflow.base_continuum
    photoionization_data = continuum.input.photoionization_data
    photoionization_index = photoionization_data.index.unique()
    upper_ion_index = pd.MultiIndex.from_arrays(
        [
            photoionization_index.get_level_values("atomic_number"),
            photoionization_index.get_level_values("ion_number") + 1,
        ],
        names=["atomic_number", "ion_number"],
    ).unique()
    radiation_field = DilutePlanckianRadiationField(
        np.asarray(plasma.t_rad) * u.K,
        np.asarray(plasma.w),
    )
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
        # The standard and IIP implementations use different cgs constant
        # sources for the complete transition table.
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
        # Independent quadrature and cgs-constant paths accumulate their
        # largest difference in the free-bound cooling channel.
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
        continuum_comparison_state.continuum.cooling_rates, process_name
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


@pytest.mark.xfail(
    strict=True,
    reason=(
        "The standalone bound-bound solver does not include IIP's coupled "
        "continuum rates, ion-stage ratios, or population-dependent Sobolev "
        "escape probabilities."
    ),
)
def test_standalone_level_populations_do_not_claim_iip_coupled_parity(
    iip_plasma_nlte_init: LegacyPlasmaArray,
) -> None:
    """Characterize the known standalone-equilibrium/IIP population gap."""
    atom_data = iip_plasma_nlte_init.atomic_data
    hydrogen_lines = atom_data.lines.loc[(1, 0, slice(None), slice(None))]
    standard_rate_matrix = RateMatrix(
        [(RadiativeRatesSolver(hydrogen_lines), "radiative")],
        atom_data.levels,
    )
    radiation_field = DilutePlanckianRadiationField(
        np.asarray(iip_plasma_nlte_init.t_rad) * u.K,
        np.asarray(iip_plasma_nlte_init.w),
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        np.asarray(iip_plasma_nlte_init.t_electrons) * u.K,
        np.asarray(iip_plasma_nlte_init.electron_densities) / u.cm**3,
    )
    standard_matrices = standard_rate_matrix.solve(
        radiation_field,
        electron_distribution,
    )
    standard_populations = (
        LevelPopulationSolver(standard_matrices, atom_data.levels)
        .solve()
        .loc[(1, 0)]
    )
    iip_populations = iip_plasma_nlte_init.level_number_density.loc[(1, 0)]
    iip_ion_population = iip_plasma_nlte_init.ion_number_density.loc[(1, 0)]
    iip_populations = iip_populations.divide(iip_ion_population, axis=1)

    common_levels = standard_populations.index.intersection(
        iip_populations.index
    )
    np.testing.assert_allclose(
        standard_populations.loc[common_levels].to_numpy(),
        iip_populations.loc[common_levels].to_numpy(),
        rtol=1e-6,
        atol=0,
    )


# comparison of plasma after the Monte Carlo calculations have been performed
def test_iip_plasma_after_mc(
    iip_regression_path,
    iip_plasma_after_mc,
):
    tau_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_tau_sobolevs_after_mc.h5",
        key="data",
    )

    beta_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_beta_sobolevs_after_mc.h5",
        key="data",
    )
    ion_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_ion_density_after_mc.h5",
        key="data",
    )
    level_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_level_number_density_after_mc.h5",
        key="data",
    )
    transition_probabilities_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_transition_probabilities_after_mc.h5",
        key="data",
    )

    # tolerances are much worse than after init

    print(
        "after MC transition_probabilities max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.transition_probabilities,
                transition_probabilities_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_after_mc.transition_probabilities,
        transition_probabilities_ctardis,
        rtol=7e-6,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "after MC ion_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.ion_number_density,
                ion_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_after_mc.ion_number_density,
        ion_number_density_ctardis,
        rtol=7e-6,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    # Sobolev values are stored differently between codes, so comparing raw data instead
    print(
        "after MC tau_sobolevs max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.tau_sobolevs,
                tau_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_after_mc.tau_sobolevs.values,
        tau_sobolevs_ctardis.values,
        rtol=7e-6,
        atol=0,
    )

    print(
        "after MC beta_sobolev max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.beta_sobolev,
                beta_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        iip_plasma_after_mc.beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=7e-6,
        atol=0,
    )

    print(
        "after MC level_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                iip_plasma_after_mc.level_number_density,
                level_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        iip_plasma_after_mc.level_number_density,
        level_number_density_ctardis,
        rtol=7e-6,
        atol=0,
        check_dtype=False,
        check_names=False,
    )


def test_thermal_balance_solver(
    iip_regression_path, type_iip_workflow, iip_plasma_after_mc
):
    type_iip_workflow.plasma_solver = iip_plasma_after_mc
    type_iip_workflow.solve_thermal_balance()

    tau_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_tau_sobolevs_after_tb.h5",
        key="data",
    )

    beta_sobolevs_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_beta_sobolevs_after_tb.h5",
        key="data",
    )
    ion_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_ion_density_after_tb.h5",
        key="data",
    )
    level_number_density_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_level_number_density_after_tb.h5",
        key="data",
    )
    transition_probabilities_ctardis = pd.read_hdf(
        iip_regression_path / "ctardis_transition_probabilities_after_tb.h5",
        key="data",
    )

    print(
        "after thermal balance transition_probabilities max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.transition_probabilities,
                transition_probabilities_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        type_iip_workflow.plasma_solver.transition_probabilities,
        transition_probabilities_ctardis,
        rtol=7e-7,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    print(
        "after thermal balance ion_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.ion_number_density,
                ion_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        type_iip_workflow.plasma_solver.ion_number_density,
        ion_number_density_ctardis,
        rtol=6e-7,
        atol=0,
        check_dtype=False,
        check_names=False,
    )

    # Sobolev values are stored differently between codes, so comparing raw data instead
    print(
        "after thermal balance tau_sobolevs max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.tau_sobolevs,
                tau_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        type_iip_workflow.plasma_solver.tau_sobolevs.values,
        tau_sobolevs_ctardis.values,
        rtol=7e-7,
        atol=0,
    )

    print(
        "after thermal balance beta_sobolev max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.beta_sobolev,
                beta_sobolevs_ctardis,
            )
        )
    )
    np.testing.assert_allclose(
        type_iip_workflow.plasma_solver.beta_sobolev.values,
        beta_sobolevs_ctardis.values,
        rtol=7e-7,
        atol=0,
    )

    print(
        "after thermal balance level_number_density max rel diff: {:.3e}".format(
            _max_rel_diff(
                type_iip_workflow.plasma_solver.level_number_density,
                level_number_density_ctardis,
            )
        )
    )
    pd.testing.assert_frame_equal(
        type_iip_workflow.plasma_solver.level_number_density,
        level_number_density_ctardis,
        rtol=6e-7,
        atol=0,
        check_dtype=False,
        check_names=False,
    )


@pytest.mark.xfail  # JOSH: This test fails because I disabled diagonalize_ma() in the BaseContinuum object to handle multi-element sims
def test_solve_continuum_state_after_nlte_init(
    iip_regression_path, type_iip_workflow, iip_plasma_nlte_init
):
    type_iip_workflow.plasma_solver = iip_plasma_nlte_init
    type_iip_workflow.solve_continuum_state(None)

    recombination_probabilities_ctardis = pd.read_hdf(
        iip_regression_path / "continuum_recomb_transition_prob_init_nlte.h5",
        key="data",
    )

    pd.testing.assert_frame_equal(
        type_iip_workflow.base_continuum.recombination_transition_probabilities.dataframe,
        recombination_probabilities_ctardis,
        rtol=4e-3,
    )


def test_solve_montecarlo(type_iip_workflow, regression_data):
    opacity_states = type_iip_workflow.solve_opacity()
    type_iip_workflow.solve_montecarlo(opacity_states, 1000)

    type_iip_workflow.initialize_spectrum_solver()
    real_packets = type_iip_workflow.spectrum_solver.spectrum_real_packets
    expected_lum_dens = regression_data.sync_ndarray(
        real_packets.luminosity_density_lambda.value
    )

    np.testing.assert_allclose(
        real_packets.luminosity_density_lambda.value,
        expected_lum_dens,
        atol=0,
        rtol=1e-12,
    )
