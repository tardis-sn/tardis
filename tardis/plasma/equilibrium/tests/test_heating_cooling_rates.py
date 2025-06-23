from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    AdiabaticThermalRates,
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)

# Assume the following fixtures are available elsewhere:
# - photoionization_cross_sections
# - collisional_cross_sections
# - thermal_electron_distribution
# - radiation_field
# - level_population
# - ion_population
# - saha_factor
# - bound_free_heating_estimator
# - stimulated_recombination_estimator
# - heating_estimator
# - collisional_ionization_rate_coefficient
# - collisional_deexcitation_rate_coefficient
# - collisional_excitation_rate_coefficient
# - time


def test_bound_free_thermal_rates_solve(
    photoionization_cross_sections,
    level_population,
    ion_population,
    thermal_electron_distribution,
    radiation_field,
    saha_factor,
):
    rates = BoundFreeThermalRates(photoionization_cross_sections)
    heating, cooling = rates.solve(
        level_population,
        ion_population,
        thermal_electron_distribution,
        saha_factor,
        radiation_field,
    )
    assert heating is not None
    assert cooling is not None


def test_bound_free_thermal_rates_solve_with_estimators(
    photoionization_cross_sections,
    level_population,
    ion_population,
    thermal_electron_distribution,
    saha_factor,
    bound_free_heating_estimator,
    stimulated_recombination_estimator,
):
    rates = BoundFreeThermalRates(photoionization_cross_sections)
    heating, cooling = rates.solve(
        level_population,
        ion_population,
        thermal_electron_distribution,
        saha_factor,
        bound_free_heating_estimator=bound_free_heating_estimator,
        stimulated_recombination_estimator=stimulated_recombination_estimator,
    )
    assert heating is not None
    assert cooling is not None


def test_free_free_thermal_rates_heating_factor(
    ion_population,
    thermal_electron_distribution,
):
    rates = FreeFreeThermalRates()
    factor = rates.heating_factor(
        ion_population, thermal_electron_distribution.number_density
    )
    assert factor is not None


def test_free_free_thermal_rates_solve(
    heating_estimator,
    thermal_electron_distribution,
    ion_population,
):
    rates = FreeFreeThermalRates()
    heating, cooling = rates.solve(
        heating_estimator,
        thermal_electron_distribution,
        ion_population,
    )
    assert heating is not None
    assert cooling is not None


def test_collisional_ionization_thermal_rates_solve(
    photoionization_cross_sections,
    thermal_electron_distribution,
    ion_population,
    level_population,
    collisional_ionization_rate_coefficient,
    saha_factor,
):
    rates = CollisionalIonizationThermalRates(photoionization_cross_sections)
    heating, cooling = rates.solve(
        thermal_electron_distribution.number_density,
        ion_population,
        level_population,
        collisional_ionization_rate_coefficient,
        saha_factor,
    )
    assert heating is not None
    assert cooling is not None


def test_collisional_bound_thermal_rates_solve(
    collisional_cross_sections,
    thermal_electron_distribution,
    collisional_deexcitation_rate_coefficient,
    collisional_excitation_rate_coefficient,
    level_population,
):
    rates = CollisionalBoundThermalRates(collisional_cross_sections)
    heating, cooling = rates.solve(
        thermal_electron_distribution.number_density,
        collisional_deexcitation_rate_coefficient,
        collisional_excitation_rate_coefficient,
        level_population,
    )
    assert heating is not None
    assert cooling is not None


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
