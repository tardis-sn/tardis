import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import (
    assert_almost_equal,
)

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    AdiabaticThermalRates,
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.equilibrium.thermal_balance import ThermalBalanceSolver
from tardis.plasma.radiation_field import DilutePlanckianRadiationField


@pytest.fixture
def thermal_electron_distribution():
    return ThermalElectronEnergyDistribution(
        0 * u.erg,
        np.ones(24) * 9992.2722969523056 * u.K,
        np.ones(24) * 2.20676447e09 * u.cm**-3,
    )


@pytest.fixture
def radiation_field():
    return DilutePlanckianRadiationField(
        np.ones(24) * 10000 * u.K, np.ones(24) * 0.5, geometry=None
    )


@pytest.fixture()
def level_population(regression_data):
    df = pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_level_number_density.csv",
        index_col=(0, 1, 2),
    )

    df.columns = [int(col) if col.isdigit() else col for col in df.columns]
    return df


@pytest.fixture()
def ion_population(regression_data):
    df = pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_ion_number_density.csv",
        index_col=(0, 1),
    )

    df.columns = [int(col) if col.isdigit() else col for col in df.columns]
    return df


@pytest.fixture()
def level_population_ratio(regression_data):
    df = pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_level_pop_ratio.csv",
        index_col=(0, 1, 2),
    )

    df.columns = [int(col) if col.isdigit() else col for col in df.columns]
    return df


@pytest.fixture()
def bound_free_heating_estimator(regression_data):
    df = pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_bf_heating_est.csv",
        index_col=(0, 1, 2),
    )

    df.columns = [int(col) if col.isdigit() else col for col in df.columns]
    return df


@pytest.fixture()
def stimulated_recombination_estimator(regression_data):
    df = pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_stim_cooling_est.csv",
        index_col=(0, 1, 2),
    )

    df.columns = [int(col) if col.isdigit() else col for col in df.columns]
    return df


@pytest.fixture()
def collisional_ionization_rate_coefficient(regression_data):
    df = pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_coll_ion_rate_coeff.csv",
        index_col=(0, 1, 2),
    )

    df.columns = [int(col) if col.isdigit() else col for col in df.columns]
    return df


@pytest.fixture()
def collisional_excitation_rate_coefficient(regression_data):
    df = pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_coll_exc_coeff.csv",
        index_col=(0, 1, 2, 3),
    )

    df.columns = [int(col) if col.isdigit() else col for col in df.columns]
    return df


@pytest.fixture()
def collisional_deexcitation_rate_coefficient(regression_data):
    df = pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_coll_deexc_coeff.csv",
        index_col=(0, 1, 2, 3),
    )

    df.columns = [int(col) if col.isdigit() else col for col in df.columns]
    return df


@pytest.fixture()
def ctardis_lines(regression_data):
    return pd.read_csv(
        regression_data.regression_data_path
        / "testdata"
        / "thermal_data"
        / "thermal_lines.csv",
        index_col=0,
    )


def test_bound_free_thermal_rates_solve(
    nlte_atom_data,
    level_population,
    ion_population,
    thermal_electron_distribution,
    radiation_field,
    level_population_ratio,
    regression_data,
):
    rates = BoundFreeThermalRates(
        nlte_atom_data.photoionization_data.query("atomic_number == 1")
    )
    actual_heating, actual_cooling = rates.solve(
        level_population,
        ion_population,
        thermal_electron_distribution,
        level_population_ratio,
        radiation_field,
    )

    assert actual_heating is not None
    assert actual_cooling is not None
    assert isinstance(actual_heating, pd.Series)
    assert isinstance(actual_cooling, pd.Series)

    expected_heating = regression_data.sync_dataframe(
        actual_heating, key="bound_free_heating"
    )
    expected_cooling = regression_data.sync_dataframe(
        actual_cooling, key="bound_free_cooling"
    )

    pdt.assert_series_equal(actual_heating, expected_heating)
    pdt.assert_series_equal(actual_cooling, expected_cooling)


@pytest.mark.parametrize(
    "heating_rate, cooling_rate",
    [(1.2809489753862688e-06, 1.235774260367592e-06)],
)
def test_bound_free_thermal_rates_solve_with_estimators(
    nlte_atom_data,
    level_population,
    ion_population,
    thermal_electron_distribution,
    level_population_ratio,
    bound_free_heating_estimator,
    stimulated_recombination_estimator,
    heating_rate,
    cooling_rate,
    regression_data,
):
    rates = BoundFreeThermalRates(
        nlte_atom_data.photoionization_data.query("atomic_number == 1")
    )
    actual_heating, actual_cooling = rates.solve(
        level_population,
        ion_population,
        thermal_electron_distribution,
        level_population_ratio,
        bound_free_heating_estimator=bound_free_heating_estimator,
        stimulated_recombination_estimator=stimulated_recombination_estimator,
    )

    # Original parametrized assertions
    assert_almost_equal(actual_heating[0], heating_rate, decimal=10)
    assert_almost_equal(actual_cooling[0], cooling_rate, decimal=10)

    # Regression data comparison
    expected_heating = regression_data.sync_dataframe(
        actual_heating, key="bound_free_heating_with_estimators"
    )
    expected_cooling = regression_data.sync_dataframe(
        actual_cooling, key="bound_free_cooling_with_estimators"
    )

    pdt.assert_series_equal(actual_heating, expected_heating)
    pdt.assert_series_equal(actual_cooling, expected_cooling)


def test_free_free_thermal_rates_heating_factor(
    ion_population,
    thermal_electron_distribution,
    regression_data,
):
    rates = FreeFreeThermalRates()
    actual_factor = rates.heating_factor(
        ion_population,
        thermal_electron_distribution.number_density.cgs.value,
    )
    expected_factor = 4.869809426191116e18

    # Original parametrized assertion
    assert_almost_equal(actual_factor[0], expected_factor, decimal=10)

    # Regression data comparison
    expected_factor_series = regression_data.sync_dataframe(
        actual_factor, key="free_free_heating_factor"
    )

    pdt.assert_series_equal(actual_factor, expected_factor_series)


@pytest.mark.parametrize(
    "ff_heating_estimator, expected_heating_rate, expected_cooling_rate",
    [(4.89135279e-24, 2.3829164962085199e-07, 6.941664530316456e-07)],
)
def test_free_free_thermal_rates_solve(
    thermal_electron_distribution,
    ion_population,
    ff_heating_estimator,
    expected_heating_rate,
    expected_cooling_rate,
    regression_data,
):
    rates = FreeFreeThermalRates()
    actual_heating_rate, actual_cooling_rate = rates.solve(
        ff_heating_estimator,
        thermal_electron_distribution,
        ion_population,
    )

    # Original parametrized assertions
    assert_almost_equal(
        actual_heating_rate[0], expected_heating_rate, decimal=10
    )
    assert_almost_equal(
        actual_cooling_rate[0], expected_cooling_rate, decimal=10
    )

    # Regression data comparison
    expected_heating_series = regression_data.sync_dataframe(
        actual_heating_rate, key="free_free_heating_rate"
    )
    expected_cooling_series = regression_data.sync_dataframe(
        actual_cooling_rate, key="free_free_cooling_rate"
    )

    pdt.assert_series_equal(actual_heating_rate, expected_heating_series)
    pdt.assert_series_equal(actual_cooling_rate, expected_cooling_series)


@pytest.mark.parametrize(
    "heating_rate, cooling_rate",
    [(1.5946196993219911e-07, 1.6125333965984259e-07)],
)
def test_collisional_ionization_thermal_rates_solve(
    nlte_atom_data,
    thermal_electron_distribution,
    ion_population,
    level_population,
    collisional_ionization_rate_coefficient,
    level_population_ratio,
    heating_rate,
    cooling_rate,
    regression_data,
):
    rates = CollisionalIonizationThermalRates(
        nlte_atom_data.photoionization_data
    )
    actual_heating, actual_cooling = rates.solve(
        thermal_electron_distribution.number_density,
        ion_population,
        level_population,
        collisional_ionization_rate_coefficient,
        level_population_ratio,
    )

    # Original parametrized assertions
    assert_almost_equal(actual_heating[0], heating_rate)
    assert_almost_equal(actual_cooling[0], cooling_rate)

    # Regression data comparison
    expected_heating = regression_data.sync_dataframe(
        actual_heating, key="collisional_ionization_heating"
    )
    expected_cooling = regression_data.sync_dataframe(
        actual_cooling, key="collisional_ionization_cooling"
    )

    pdt.assert_series_equal(actual_heating, expected_heating)
    pdt.assert_series_equal(actual_cooling, expected_cooling)


@pytest.mark.parametrize(
    "heating_rate, cooling_rate",
    [(8.49837495539e-06, 8.5059159013e-06)],
)
def test_collisional_bound_thermal_rates_solve(
    ctardis_lines,
    thermal_electron_distribution,
    collisional_deexcitation_rate_coefficient,
    collisional_excitation_rate_coefficient,
    level_population,
    heating_rate,
    cooling_rate,
    regression_data,
):
    rates = CollisionalBoundThermalRates(ctardis_lines)
    actual_heating, actual_cooling = rates.solve(
        thermal_electron_distribution.number_density,
        collisional_deexcitation_rate_coefficient,
        collisional_excitation_rate_coefficient,
        level_population,
    )

    # Original parametrized assertions
    assert_almost_equal(actual_heating[0], heating_rate, decimal=8)
    assert_almost_equal(actual_cooling[0], cooling_rate, decimal=8)

    # Regression data comparison
    expected_heating = regression_data.sync_dataframe(
        actual_heating, key="collisional_bound_heating"
    )
    expected_cooling = regression_data.sync_dataframe(
        actual_cooling, key="collisional_bound_cooling"
    )

    pdt.assert_series_equal(actual_heating, expected_heating)
    pdt.assert_series_equal(actual_cooling, expected_cooling)


@pytest.mark.parametrize(
    "time, expected_cooling_rate",
    [(1 * u.s, 0.009133236799630136 * u.erg / (u.s * u.cm**3))],
)
def test_adiabatic_thermal_rates_solve(
    thermal_electron_distribution,
    time,
    expected_cooling_rate,
):
    rates = AdiabaticThermalRates()
    actual_cooling_rate = rates.solve(
        thermal_electron_distribution,
        time,
    )

    # Original parametrized assertion
    assert_quantity_allclose(
        actual_cooling_rate, expected_cooling_rate, rtol=1e-14
    )


def test_thermal_balance_solver(
    thermal_electron_distribution,
    radiation_field,
    level_population,
    ion_population,
    collisional_ionization_rate_coefficient,
    collisional_deexcitation_rate_coefficient,
    collisional_excitation_rate_coefficient,
    stimulated_recombination_estimator,
    bound_free_heating_estimator,
    level_population_ratio,
    ctardis_lines,
    nlte_atom_data,
    regression_data,
):
    bound_free_rates_solver = BoundFreeThermalRates(
        nlte_atom_data.photoionization_data.query("atomic_number == 1")
    )
    free_free_rates_solver = FreeFreeThermalRates()
    collisional_ionization_rates_solver = CollisionalIonizationThermalRates(
        nlte_atom_data.photoionization_data.query("atomic_number == 1")
    )
    collisional_bound_rates_solver = CollisionalBoundThermalRates(ctardis_lines)

    ff_heating_estimator = 4.89135279e-24  # from chvogl's code

    thermal_balance_solver = ThermalBalanceSolver(
        bound_free_rates_solver,
        free_free_rates_solver,
        collisional_ionization_rates_solver,
        collisional_bound_rates_solver,
    )

    actual_total_heating_rate, actual_fractional_heating_rate = (
        thermal_balance_solver.solve(
            thermal_electron_distribution,
            level_population,
            ion_population,
            collisional_ionization_rate_coefficient,
            collisional_deexcitation_rate_coefficient,
            collisional_excitation_rate_coefficient,
            ff_heating_estimator,
            level_population_ratio,
            radiation_field,
            bound_free_heating_estimator,
            stimulated_recombination_estimator,
        )
    )

    expected_total_heating_rate = -4.2003240358632565e-07
    expected_fractional_heating_rate = -0.03963650747926685

    # Original parametrized assertions
    assert_almost_equal(
        actual_total_heating_rate[0], expected_total_heating_rate, decimal=14
    )
    assert_almost_equal(
        actual_fractional_heating_rate[0],
        expected_fractional_heating_rate,
        decimal=14,
    )

    # Regression data comparison
    expected_total_heating_series = regression_data.sync_dataframe(
        actual_total_heating_rate, key="thermal_balance_total_heating_rate"
    )
    expected_fractional_heating_series = regression_data.sync_dataframe(
        actual_fractional_heating_rate,
        key="thermal_balance_fractional_heating_rate",
    )

    pdt.assert_series_equal(
        actual_total_heating_rate, expected_total_heating_series
    )
    pdt.assert_series_equal(
        actual_fractional_heating_rate, expected_fractional_heating_series
    )
