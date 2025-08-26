import astropy.units as u
import numpy as np
import pandas as pd
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
from tardis.plasma.radiation_field import DilutePlanckianRadiationField


@pytest.fixture
def thermal_electron_distribution():
    return ThermalElectronEnergyDistribution(
        0 * u.erg, 9992.2722969523056 * u.K, 2.20676447e09 * u.cm**-3
    )


@pytest.fixture
def radiation_field():
    return DilutePlanckianRadiationField(
        np.ones(1) * 10000 * u.K, np.array([0.5]), geometry=None
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
):
    rates = BoundFreeThermalRates(
        nlte_atom_data.photoionization_data.query("atomic_number == 1")
    )
    heating, cooling = rates.solve(
        level_population[0],
        ion_population[0],
        thermal_electron_distribution,
        level_population_ratio[0],
        radiation_field,
    )
    assert heating is not None
    assert cooling is not None


@pytest.mark.parametrize(
    "heating_rate, cooling_rate",
    [(1.2809489753862688e-06, 1.2018593543520837e-06)],
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
):
    rates = BoundFreeThermalRates(
        nlte_atom_data.photoionization_data.query("atomic_number == 1")
    )
    heating, cooling = rates.solve(
        level_population[0],
        ion_population[0],
        thermal_electron_distribution,
        level_population_ratio[0],
        bound_free_heating_estimator=bound_free_heating_estimator[0],
        stimulated_recombination_estimator=stimulated_recombination_estimator[
            0
        ],
    )

    assert_almost_equal(heating, heating_rate)
    assert_almost_equal(cooling, cooling_rate)


def test_free_free_thermal_rates_heating_factor(
    ion_population,
    thermal_electron_distribution,
):
    rates = FreeFreeThermalRates()
    factor = rates.heating_factor(
        ion_population[0],
        thermal_electron_distribution.number_density.cgs.value,
    )
    expected_factor = 4.869809426186787e18

    assert_almost_equal(factor, expected_factor)


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
):
    rates = FreeFreeThermalRates()
    actual_heating_rate, actual_cooling_rate = rates.solve(
        ff_heating_estimator,
        thermal_electron_distribution,
        ion_population[0],  # only one cell
    )
    assert_almost_equal(actual_heating_rate, expected_heating_rate)
    assert_almost_equal(actual_cooling_rate, expected_cooling_rate)


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
):
    rates = CollisionalIonizationThermalRates(
        nlte_atom_data.photoionization_data
    )
    heating, cooling = rates.solve(
        thermal_electron_distribution.number_density,
        ion_population[0],
        level_population[0],
        collisional_ionization_rate_coefficient[0],
        level_population_ratio[0],
    )
    assert_almost_equal(heating, heating_rate)
    assert_almost_equal(cooling, cooling_rate)


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
):
    rates = CollisionalBoundThermalRates(ctardis_lines)
    heating, cooling = rates.solve(
        thermal_electron_distribution.number_density,
        collisional_deexcitation_rate_coefficient[0],
        collisional_excitation_rate_coefficient[0],
        level_population[0],
    )
    assert_almost_equal(heating, heating_rate)
    assert_almost_equal(cooling, cooling_rate)


@pytest.mark.parametrize(
    "time, expected_cooling_rate",
    [(1 * u.s, 0.009133236799630136 * u.erg / (u.s * u.cm**3))],
)
def test_adiabatic_thermal_rates_solve(
    thermal_electron_distribution, time, expected_cooling_rate
):
    rates = AdiabaticThermalRates()
    actual_cooling_rate = rates.solve(
        thermal_electron_distribution,
        time,
    )
    assert_quantity_allclose(
        actual_cooling_rate, expected_cooling_rate, rtol=1e-14
    )
