from pathlib import Path

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis import constants as const
from tardis.iip_plasma.properties.continuum import (
    BfHeatingRateCoeff,
    PhotoIonRateCoeff,
    SpontRecombRateCoeff,
    StimRecombRateCoeff,
    ThermalBalanceTest,
)
from tardis.io.atom_data import AtomData
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
)
from tardis.plasma.equilibrium.rates.photoionization_rates import (
    EstimatedPhotoionizationRateSolver,
)
from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_ionization_rate_dataframe,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.transport.montecarlo.estimators import init_estimators_continuum

REFERENCE_RADIATION_TEMPERATURES_K = np.array([10000.0, 12000.0])
REFERENCE_RADIATION_DILUTION_FACTORS = np.array([0.4, 0.8])
REFERENCE_ELECTRON_TEMPERATURES_K = np.array([9000.0, 11000.0])
REFERENCE_ELECTRON_DENSITY_CM3 = 1.0e9
REFERENCE_ESTIMATOR_TIME_S = 2.0e5
REFERENCE_ESTIMATOR_VOLUME_CM3 = 3.0e30


@pytest.fixture
def photoionization_data(nlte_atom_data: AtomData) -> pd.DataFrame:
    return nlte_atom_data.photoionization_data.loc[(1, 0, [1])].sort_values(
        ["atomic_number", "ion_number", "level_number", "nu"]
    )


@pytest.fixture
def lyman_photoionization_data(nlte_atom_data: AtomData) -> pd.DataFrame:
    return nlte_atom_data.photoionization_data.loc[(1, 0, [0, 1])].sort_values(
        ["atomic_number", "ion_number", "level_number", "nu"]
    )


@pytest.fixture
def radiation_field() -> DilutePlanckianRadiationField:
    return DilutePlanckianRadiationField(
        REFERENCE_RADIATION_TEMPERATURES_K * u.K,
        REFERENCE_RADIATION_DILUTION_FACTORS,
    )


@pytest.fixture
def analytic_photoionization_rates(
    photoionization_data: pd.DataFrame,
    radiation_field: DilutePlanckianRadiationField,
) -> dict[str, pd.DataFrame]:
    photo_data = photoionization_data
    electron_temperature = REFERENCE_ELECTRON_TEMPERATURES_K * u.K
    standard_solver = AnalyticPhotoionizationCoeffSolver(photo_data)
    gamma, alpha_stim = standard_solver.solve(
        radiation_field, electron_temperature
    )

    iip_j_nu = PhotoIonRateCoeff._calculate_j_nus(
        photo_data,
        radiation_field.dilution_factor,
        radiation_field.temperature_kelvin,
    )
    standard_j_nu = standard_solver.calculate_mean_intensity_photoionization_df(
        radiation_field
    )
    iip_gamma = PhotoIonRateCoeff(None).calculate_from_radiation_field_model(
        photo_data,
        radiation_field.dilution_factor,
        radiation_field.temperature_kelvin,
    )
    iip_alpha_stim = StimRecombRateCoeff(
        None
    ).calculate_from_radiation_field_model(
        photo_data,
        radiation_field.dilution_factor,
        radiation_field.temperature_kelvin,
        None,
        electron_temperature.value,
        pd.DataFrame(
            np.ones((1, 2)),
            index=photo_data.index.unique(),
            columns=[0, 1],
        ),
    )
    return {
        "standard_j_nu": standard_j_nu,
        "iip_j_nu": iip_j_nu,
        "gamma": gamma,
        "iip_gamma": iip_gamma,
        "alpha_stim": alpha_stim,
        "iip_alpha_stim": iip_alpha_stim,
    }


def test_photoionization_mean_intensity_matches_iip(
    analytic_photoionization_rates: dict[str, pd.DataFrame],
) -> None:
    npt.assert_allclose(
        analytic_photoionization_rates["standard_j_nu"].to_numpy(),
        analytic_photoionization_rates["iip_j_nu"].to_numpy(),
        rtol=1e-12,
    )


def test_photoionization_rate_matches_iip(
    analytic_photoionization_rates: dict[str, pd.DataFrame],
) -> None:
    npt.assert_allclose(
        analytic_photoionization_rates["gamma"].to_numpy(),
        analytic_photoionization_rates["iip_gamma"].to_numpy(),
        rtol=2e-7,
    )


def test_stimulated_recombination_rate_matches_iip(
    analytic_photoionization_rates: dict[str, pd.DataFrame],
) -> None:
    npt.assert_allclose(
        analytic_photoionization_rates["alpha_stim"].to_numpy(),
        analytic_photoionization_rates["iip_alpha_stim"].to_numpy(),
        rtol=2e-6,
    )


def test_zero_radiation_gives_zero_photoionization_and_stimulated_recombination(
    photoionization_data: pd.DataFrame,
) -> None:
    photo_data = photoionization_data
    zero_field = DilutePlanckianRadiationField(
        REFERENCE_RADIATION_TEMPERATURES_K * u.K, np.zeros(2)
    )
    gamma, alpha_stim = AnalyticPhotoionizationCoeffSolver(photo_data).solve(
        zero_field, REFERENCE_ELECTRON_TEMPERATURES_K * u.K
    )
    assert np.all(gamma.to_numpy() == 0.0)
    assert np.all(alpha_stim.to_numpy() == 0.0)


def test_spontaneous_recombination_is_positive_and_lyman_suppression_is_explicit(
    lyman_photoionization_data: pd.DataFrame,
) -> None:
    photo_data = lyman_photoionization_data
    temperatures = REFERENCE_RADIATION_TEMPERATURES_K * u.K
    standard_alpha = SpontaneousRecombinationCoeffSolver(photo_data).solve(
        temperatures
    )
    assert (standard_alpha.loc[(1, 0, 1)] >= 0).all()
    assert np.all(standard_alpha.loc[(1, 0, 0)] == 0.0)

    phi_lucy = pd.DataFrame(
        np.ones((2, 2)),
        index=photo_data.index.unique(),
        columns=[0, 1],
    )
    iip_alpha = SpontRecombRateCoeff(
        type("IterationState", (), {"niter": 2, "niter_ly": 1})()
    ).calculate(photo_data, temperatures.value, phi_lucy)
    assert np.all(iip_alpha.loc[(1, 0, 0)] == 0.0)
    npt.assert_allclose(
        standard_alpha.loc[(1, 0, 1)].to_numpy(),
        iip_alpha.loc[(1, 0, 1)].to_numpy(),
        rtol=2e-4,
    )


def test_estimator_coefficients_reproduce_regression_inputs(
    tardis_regression_path: Path,
) -> None:
    edge_index = pd.MultiIndex.from_tuples(
        [(1, 0, 1), (1, 0, 2)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    photo_ion_estimator = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_photo_ion_estimator_after_mc.h5",
        key="data",
    ).loc[edge_index, [0, 1]]
    stim_recomb_estimator = pd.read_hdf(
        tardis_regression_path
        / "tardis"
        / "workflows"
        / "tests"
        / "ctardis_stim_recomb_estimator_after_mc.h5",
        key="data",
    ).loc[edge_index, [0, 1]]
    estimators = init_estimators_continuum(
        n_levels_bf_species_by_n_cells_tuple=(2, 2), n_cells=2
    )
    estimators.photo_ion_estimator[:] = photo_ion_estimator.to_numpy()
    estimators.stim_recomb_estimator[:] = stim_recomb_estimator.to_numpy()
    time_simulation = REFERENCE_ESTIMATOR_TIME_S * u.s
    volume = REFERENCE_ESTIMATOR_VOLUME_CM3 * u.cm**3
    normalization = 1.0 / (
        time_simulation.to_value(u.s)
        * volume.to_value(u.cm**3)
        * const.h.cgs.value
    )
    gamma, alpha_stim = EstimatedPhotoionizationCoeffSolver(
        pd.Series([0, 1], index=edge_index)
    ).solve(estimators, time_simulation, volume)
    npt.assert_allclose(
        gamma.to_numpy(),
        photo_ion_estimator.to_numpy() * normalization,
    )
    npt.assert_allclose(
        alpha_stim.to_numpy(),
        stim_recomb_estimator.to_numpy() * normalization,
    )
    assert gamma.index.equals(edge_index)


def test_estimated_rates_use_lucy_ion_matrix_coefficients(
    mock_photoionization_cross_sections: pd.DataFrame,
) -> None:
    """Estimated rates use level fractions and Lucy's Saha factor."""
    level_index = mock_photoionization_cross_sections.index
    columns = pd.Index([0])
    level_population = pd.DataFrame(
        [2.0, 3.0], index=level_index, columns=columns
    )
    ion_population = pd.DataFrame(
        [10.0, 20.0],
        index=pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
        ),
        columns=columns,
    )
    level_to_continuum_saha_factor = pd.DataFrame(
        [0.5, 0.25], index=level_index, columns=columns
    )
    estimators = init_estimators_continuum(
        n_levels_bf_species_by_n_cells_tuple=(2, 1), n_cells=1
    )
    estimators.photo_ion_estimator[:] = [[2.0], [4.0]]
    estimators.stim_recomb_estimator[:] = [[0.1], [0.2]]
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg, np.array([1.0e4]) * u.K, np.array([4.0]) / u.cm**3
    )

    solver = EstimatedPhotoionizationRateSolver(
        mock_photoionization_cross_sections,
        pd.Series([0, 1], index=level_index),
        estimators,
        1 * u.s,
        1 * u.cm**3,
    )
    actual_photoionization, actual_recombination = solver.solve(
        electron_distribution,
        level_population,
        ion_population,
        level_to_continuum_saha_factor,
    )

    coefficient_solver = EstimatedPhotoionizationCoeffSolver(
        pd.Series([0, 1], index=level_index)
    )
    photoionization_coefficient, stimulated_coefficient = (
        coefficient_solver.solve(estimators, 1 * u.s, 1 * u.cm**3)
    )
    spontaneous_coefficient = SpontaneousRecombinationCoeffSolver(
        mock_photoionization_cross_sections
    ).solve(electron_distribution.temperature)
    expected_photoionization = pd.DataFrame(
        np.asarray(photoionization_coefficient),
        index=level_index,
        columns=columns,
    ) * pd.DataFrame([0.2, 0.3], index=level_index, columns=columns)
    expected_photoionization.loc[(1, 0, 0)] = 0.0
    spontaneous_recomb_coeff_df = pd.DataFrame(
        np.asarray(spontaneous_coefficient),
        index=level_index,
        columns=columns,
    )
    stimulated_recomb_coeff_df = pd.DataFrame(
        np.asarray(stimulated_coefficient),
        index=level_index,
        columns=columns,
    )
    expected_recombination = (
        (spontaneous_recomb_coeff_df + stimulated_recomb_coeff_df)
        * level_to_continuum_saha_factor
        * 4.0
    )
    expected_recombination.loc[(1, 0, 0)] = 0.0
    expected_photoionization = reindex_ionization_rate_dataframe(
        expected_photoionization, recombination=False
    )
    expected_recombination = reindex_ionization_rate_dataframe(
        expected_recombination, recombination=True
    )

    assert actual_photoionization.index.equals(expected_photoionization.index)
    assert actual_recombination.index.equals(expected_recombination.index)
    npt.assert_allclose(
        actual_photoionization.to_numpy(), expected_photoionization.to_numpy()
    )
    npt.assert_allclose(
        actual_recombination.to_numpy(), expected_recombination.to_numpy()
    )


def test_bound_free_heating_and_cooling_match_iip_plasma(
    lyman_photoionization_data: pd.DataFrame,
) -> None:
    photo_data = lyman_photoionization_data
    temperatures = REFERENCE_RADIATION_TEMPERATURES_K
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        temperatures * u.K,
        np.ones(2) * REFERENCE_ELECTRON_DENSITY_CM3 / u.cm**3,
    )
    level_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    level_population = pd.DataFrame(
        np.ones((2, 2)), index=level_index, columns=[0, 1]
    )
    ion_population = pd.DataFrame(
        np.ones((1, 2)),
        index=pd.MultiIndex.from_tuples(
            [(1, 1)], names=["atomic_number", "ion_number"]
        ),
        columns=[0, 1],
    )
    level_population_ratio = pd.DataFrame(
        np.ones((2, 2)), index=level_index, columns=[0, 1]
    )
    # The IIP workflow supplies bound-free heating through its MC estimator.
    # Its Lyman-continuum heating coefficient is suppressed in this state.
    bf_heating_estimator = pd.DataFrame(
        [[0.0, 0.0], [3.0, 4.0]], index=level_index, columns=[0, 1]
    )
    heating, cooling = BoundFreeThermalRates(photo_data).solve(
        level_population,
        ion_population,
        electron_distribution,
        level_population_ratio,
        bound_free_heating_estimator=bf_heating_estimator,
    )
    iip_heating_coeff = BfHeatingRateCoeff(
        type("IterationState", (), {"niter": 0, "niter_ly": -1})()
    ).calculate(
        bf_heating_estimator,
        level_index,
    )
    iip_heating = (
        iip_heating_coeff * level_population.loc[iip_heating_coeff.index]
    ).sum()

    thermal_balance = ThermalBalanceTest(None)
    iip_cooling_coeff = pd.DataFrame(
        {
            cell: thermal_balance._calculate_sp_recomb_heating_rate_coeff(
                temperature, photo_data
            )
            for cell, temperature in enumerate(temperatures)
        }
    )
    iip_cooling_coeff.loc[(1, 0, 0)] = 0.0
    iip_cooling = (
        iip_cooling_coeff
        * level_population_ratio.loc[iip_cooling_coeff.index]
        * electron_distribution.number_density.value
        * ion_population.loc[(1, 1)]
    ).sum()

    npt.assert_allclose(heating.to_numpy(), iip_heating.to_numpy(), rtol=1e-12)
    npt.assert_allclose(cooling.to_numpy(), iip_cooling.to_numpy(), rtol=2e-6)


def test_bound_free_non_estimator_rates_match_iip_plasma(
    lyman_photoionization_data: pd.DataFrame,
    radiation_field: DilutePlanckianRadiationField,
) -> None:
    photo_data = lyman_photoionization_data
    temperatures = REFERENCE_RADIATION_TEMPERATURES_K
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg,
        temperatures * u.K,
        np.ones(2) * REFERENCE_ELECTRON_DENSITY_CM3 / u.cm**3,
    )
    level_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    level_population = pd.DataFrame(
        np.ones((2, 2)), index=level_index, columns=[0, 1]
    )
    ion_population = pd.DataFrame(
        np.ones((1, 2)),
        index=pd.MultiIndex.from_tuples(
            [(1, 1)], names=["atomic_number", "ion_number"]
        ),
        columns=[0, 1],
    )
    level_population_ratio = pd.DataFrame(
        np.ones((2, 2)), index=level_index, columns=[0, 1]
    )
    heating, cooling = BoundFreeThermalRates(photo_data).solve(
        level_population,
        ion_population,
        electron_distribution,
        level_population_ratio,
        radiation_field,
    )

    iip_heating_coeff = BfHeatingRateCoeff(
        None
    ).calculate_from_radiation_field_model(
        photo_data,
        radiation_field.dilution_factor,
        radiation_field.temperature_kelvin,
    )
    iip_heating = (
        iip_heating_coeff * level_population.loc[iip_heating_coeff.index]
    ).sum()

    thermal_balance = ThermalBalanceTest(None)
    iip_cooling_coeff = pd.DataFrame(
        {
            cell: thermal_balance._calculate_sp_recomb_heating_rate_coeff(
                temperature, photo_data
            )
            for cell, temperature in enumerate(temperatures)
        }
    )
    iip_cooling_coeff.loc[(1, 0, 0)] = 0.0
    iip_cooling = (
        iip_cooling_coeff
        * level_population_ratio.loc[iip_cooling_coeff.index]
        * electron_distribution.number_density.value
        * ion_population.loc[(1, 1)]
    ).sum()

    npt.assert_allclose(heating.to_numpy(), iip_heating.to_numpy(), rtol=2e-2)
    npt.assert_allclose(cooling.to_numpy(), iip_cooling.to_numpy(), rtol=2e-6)
