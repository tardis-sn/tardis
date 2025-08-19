import astropy.units as u
import numpy as np
import pytest
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
def level_population():
    return None


@pytest.fixture()
def ion_population():
    return None


@pytest.fixture()
def saha_factor():
    return None


@pytest.fixture()
def bound_free_heating_estimator():
    return None


@pytest.fixture()
def stimulated_recombination_estimator():
    return None


@pytest.fixture()
def collisional_ionization_rate_coefficient():
    return None


@pytest.fixture()
def collisional_excitation_rate_coefficient():
    return None


@pytest.fixture()
def collisional_deexcitation_rate_coefficient():
    return None


def test_bound_free_thermal_rates_solve(
    nlte_atom_data,
    level_population,
    ion_population,
    thermal_electron_distribution,
    radiation_field,
    saha_factor,
):
    rates = BoundFreeThermalRates(nlte_atom_data.photoionization_data)
    heating, cooling = rates.solve(
        level_population,
        ion_population,
        thermal_electron_distribution,
        saha_factor,
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
    saha_factor,
    bound_free_heating_estimator,
    stimulated_recombination_estimator,
    heating_rate,
    cooling_rate,
):
    rates = BoundFreeThermalRates(nlte_atom_data.photoionization_data)
    heating, cooling = rates.solve(
        level_population,
        ion_population,
        thermal_electron_distribution,
        saha_factor,
        bound_free_heating_estimator=bound_free_heating_estimator,
        stimulated_recombination_estimator=stimulated_recombination_estimator,
    )

    assert_almost_equal(heating, heating_rate)
    assert_almost_equal(cooling, cooling_rate)


def test_free_free_thermal_rates_heating_factor(
    ion_population,
    thermal_electron_distribution,
):
    rates = FreeFreeThermalRates()
    factor = rates.heating_factor(
        ion_population, thermal_electron_distribution.number_density
    )
    assert factor is not None


@pytest.mark.parametrize(
    "heating_rate, cooling_rate",
    [(2.3829164962085199e-07, 6.941664530316456e-07)],
)
def test_free_free_thermal_rates_solve(
    heating_estimator,
    thermal_electron_distribution,
    ion_population,
    heating_rate,
    cooling_rate,
):
    rates = FreeFreeThermalRates()
    heating, cooling = rates.solve(
        heating_estimator,
        thermal_electron_distribution,
        ion_population,
    )
    assert_almost_equal(heating, heating_rate)
    assert_almost_equal(cooling, cooling_rate)


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
    saha_factor,
    heating_rate,
    cooling_rate,
):
    rates = CollisionalIonizationThermalRates(
        nlte_atom_data.photoionization_data
    )
    heating, cooling = rates.solve(
        thermal_electron_distribution.number_density,
        ion_population,
        level_population,
        collisional_ionization_rate_coefficient,
        saha_factor,
    )
    assert_almost_equal(heating, heating_rate)
    assert_almost_equal(cooling, cooling_rate)


@pytest.mark.parametrize(
    "heating_rate, cooling_rate",
    [(8.49837495539e-06, 8.5059159013e-06)],
)
def test_collisional_bound_thermal_rates_solve(
    lines,
    thermal_electron_distribution,
    collisional_deexcitation_rate_coefficient,
    collisional_excitation_rate_coefficient,
    level_population,
    heating_rate,
    cooling_rate,
):
    rates = CollisionalBoundThermalRates(lines)
    heating, cooling = rates.solve(
        thermal_electron_distribution.number_density,
        collisional_deexcitation_rate_coefficient,
        collisional_excitation_rate_coefficient,
        level_population,
    )
    assert_almost_equal(heating, heating_rate)
    assert_almost_equal(cooling, cooling_rate)


def test_adiabatic_thermal_rates_solve(
    thermal_electron_distribution,
    time,
):
    rates = AdiabaticThermalRates()
    cooling = rates.solve(
        thermal_electron_distribution,
        time,
    )
    assert cooling is not None
