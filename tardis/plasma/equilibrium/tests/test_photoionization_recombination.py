from pathlib import Path

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u
from tardisbase.testing.regression_data.regression_data import RegressionData

from tardis import constants as const
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
def helium_photoionization_data() -> pd.DataFrame:
    """Two helium levels, each sampled at two continuum frequencies."""
    level_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0), (2, 0, 0), (2, 0, 1), (2, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    return pd.DataFrame(
        {
            "nu": [1.0e15, 2.0e15, 3.0e15, 4.0e15],
            "x_sect": [1.0e-18] * 4,
        },
        index=level_index,
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
    standard_j_nu = standard_solver.calculate_mean_intensity_photoionization_df(
        radiation_field
    )

    return {
        "standard_j_nu": standard_j_nu,
        "gamma": gamma,
        "alpha_stim": alpha_stim,
    }


def test_photoionization_mean_intensity_matches_iip(
    analytic_photoionization_rates: dict[str, pd.DataFrame],
    regression_data: RegressionData,
) -> None:
    actual = pd.DataFrame(
        analytic_photoionization_rates["standard_j_nu"].to_numpy()
    )
    expected = regression_data.sync_dataframe(actual, key="allclose_0")
    npt.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=1e-12)


def test_photoionization_rate_matches_iip(
    analytic_photoionization_rates: dict[str, pd.DataFrame],
    regression_data: RegressionData,
) -> None:
    actual = pd.DataFrame(analytic_photoionization_rates["gamma"].to_numpy())
    expected = regression_data.sync_dataframe(actual, key="allclose_0")
    npt.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=2e-7)


def test_stimulated_recombination_rate_matches_iip(
    analytic_photoionization_rates: dict[str, pd.DataFrame],
    regression_data: RegressionData,
) -> None:
    actual = pd.DataFrame(
        analytic_photoionization_rates["alpha_stim"].to_numpy()
    )
    expected = regression_data.sync_dataframe(actual, key="allclose_0")
    npt.assert_allclose(actual.to_numpy(), expected.to_numpy(), rtol=2e-6)


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
    regression_data: RegressionData,
) -> None:
    photo_data = lyman_photoionization_data
    temperatures = REFERENCE_RADIATION_TEMPERATURES_K * u.K
    standard_alpha = SpontaneousRecombinationCoeffSolver(photo_data).solve(
        temperatures
    )
    assert (standard_alpha.loc[(1, 0, 1)] >= 0).all()
    assert np.all(standard_alpha.loc[(1, 0, 0)] == 0.0)

    actual = pd.DataFrame(standard_alpha.loc[(1, 0, 1)].to_numpy())
    expected = regression_data.sync_dataframe(actual, key="allclose_0")
    npt.assert_allclose(
        actual.to_numpy().ravel(), expected.to_numpy().ravel(), rtol=2e-4
    )


def test_helium_spontaneous_recombination_keeps_input_levels(
    helium_photoionization_data: pd.DataFrame,
) -> None:
    """Helium-only continua produce only helium recombination coefficients.

    Claim
    -----
    Each available helium level receives its spontaneous coefficient; an
    absent hydrogen level is not introduced by Case B handling.

    Regime
    ------
    Two He I levels with constant cross section sampled at two frequencies.

    Verification
    ------------
    Compare against the explicitly simplified Milne-integrand trapezoid.
    """
    temperature = 1.0e4
    actual = SpontaneousRecombinationCoeffSolver(
        helium_photoionization_data
    ).solve(np.array([temperature]) * u.K)

    expected = []
    for _, cross_section in helium_photoionization_data.groupby(level=[0, 1, 2]):
        frequencies = cross_section.nu.to_numpy()
        integrand = (
            8.0
            * np.pi
            * cross_section.x_sect.to_numpy()
            * frequencies**2
            / const.c.cgs.value**2
            * np.exp(
                -const.h.cgs.value
                * frequencies
                / (const.k_B.cgs.value * temperature)
            )
        )
        expected.append(np.trapezoid(integrand, frequencies))

    pdt.assert_index_equal(
        actual.index, helium_photoionization_data.index.unique()
    )
    npt.assert_allclose(
        actual.iloc[:, 0].to_numpy(), expected, rtol=1.0e-12
    )


def test_estimated_helium_rates_do_not_require_hydrogen_continuum(
    helium_photoionization_data: pd.DataFrame,
) -> None:
    """The estimated rate matrix accepts helium-only continuum data.

    Claim
    -----
    Photoionization and recombination rates preserve the supplied He I levels.

    Regime
    ------
    One shell, two He I levels, and fixed packet estimators.

    Verification
    ------------
    Independently apply level fractions and the Lucy Saha/electron factors.
    """
    level_index = helium_photoionization_data.index.unique()
    columns = pd.Index([0])
    estimators = init_estimators_continuum((2, 1), 1)
    estimators.photo_ion_estimator[:] = (
        np.array([[2.0], [3.0]]) * const.h.cgs.value
    )
    estimators.stim_recomb_estimator[:] = (
        np.array([[0.1], [0.2]]) * const.h.cgs.value
    )
    level_population = pd.DataFrame(
        [[2.0], [3.0]], index=level_index, columns=columns
    )
    ion_population = pd.DataFrame(
        [[10.0], [5.0]],
        index=pd.MultiIndex.from_tuples(
            [(2, 0), (2, 1)], names=["atomic_number", "ion_number"]
        ),
        columns=columns,
    )
    saha_factor = pd.DataFrame(
        [[0.5], [0.25]], index=level_index, columns=columns
    )
    electron_distribution = ThermalElectronEnergyDistribution(
        0 * u.erg, np.array([1.0e4]) * u.K, np.array([4.0]) / u.cm**3
    )
    solver = EstimatedPhotoionizationRateSolver(
        helium_photoionization_data,
        pd.Series([0, 1], index=level_index),
        estimators,
        1.0 * u.s,
        1.0 * u.cm**3,
    )

    photoionization, recombination = solver.solve(
        electron_distribution,
        level_population,
        ion_population,
        saha_factor,
    )

    spontaneous = np.array(
        [
            np.trapezoid(
                8.0
                * np.pi
                * group.x_sect.to_numpy()
                * group.nu.to_numpy() ** 2
                / const.c.cgs.value**2
                * np.exp(
                    -const.h.cgs.value
                    * group.nu.to_numpy()
                    / (const.k_B.cgs.value * 1.0e4)
                ),
                group.nu.to_numpy(),
            )
            for _, group in helium_photoionization_data.groupby(level=[0, 1, 2])
        ]
    )
    npt.assert_allclose(photoionization.iloc[:, 0], [0.4, 0.9], rtol=1.0e-12)
    npt.assert_allclose(
        recombination.iloc[:, 0],
        (spontaneous + np.array([0.1, 0.2])) * np.array([0.5, 0.25]) * 4.0,
        rtol=1.0e-12,
    )
    assert photoionization.index.get_level_values("atomic_number").unique().tolist() == [2]
    assert recombination.index.get_level_values("atomic_number").unique().tolist() == [2]


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
    regression_data: RegressionData,
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
    expected_heating = regression_data.sync_dataframe(
        pd.DataFrame(heating.to_numpy()), key="allclose_0"
    )
    expected_cooling = regression_data.sync_dataframe(
        pd.DataFrame(cooling.to_numpy()), key="allclose_1"
    )
    npt.assert_allclose(
        heating.to_numpy(), expected_heating.to_numpy().ravel(), rtol=1e-12
    )
    npt.assert_allclose(
        cooling.to_numpy(), expected_cooling.to_numpy().ravel(), rtol=2e-6
    )


def test_bound_free_non_estimator_rates_match_iip_plasma(
    lyman_photoionization_data: pd.DataFrame,
    radiation_field: DilutePlanckianRadiationField,
    regression_data: RegressionData,
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

    expected_heating = regression_data.sync_dataframe(
        pd.DataFrame(heating.to_numpy()), key="allclose_0"
    )
    expected_cooling = regression_data.sync_dataframe(
        pd.DataFrame(cooling.to_numpy()), key="allclose_1"
    )
    npt.assert_allclose(
        heating.to_numpy(), expected_heating.to_numpy().ravel(), rtol=2e-2
    )
    npt.assert_allclose(
        cooling.to_numpy(), expected_cooling.to_numpy().ravel(), rtol=2e-6
    )
