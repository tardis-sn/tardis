import numpy as np

from tardis import constants as const
from tardis.transport.montecarlo.estimators.util import (
    integrate_array_by_blocks,
)


class BoundFreeThermalRates:
    """
    Class to represent the bound-free heating rate.
    """

    def __init__(self, photoionization_cross_sections):
        self.photoionization_cross_sections = photoionization_cross_sections
        self.nu = photoionization_cross_sections["nu"].values
        self.photoionization_block_references = np.pad(
            self.photoionization_cross_sections.nu.groupby(level=[0, 1, 2])
            .count()
            .values.cumsum(),
            [1, 0],
        )

    def solve(
        self,
        level_number_density,
        ion_number_density,
        electron_density,
        radiation_field,
        boltzmann_factor,
        phi,
        bound_free_heating_estimator=None,
        stimulated_recombination_estimator=None,
    ):
        nu_i = self.nu.groupby(level=[0, 1, 2]).first()
        nu_is = nu_i.loc[self.photoionization_cross_sections.index]
        mean_intensities = radiation_field.calculate_mean_intensity(self.nu)

        if bound_free_heating_estimator is not None:
            # TODO: check if this is correct
            integrated_heating_coefficient = bound_free_heating_estimator
        else:
            heating_coefficient = (
                (
                    4
                    * np.pi
                    * self.photoionization_cross_sections["x_sect"]
                    * self.nu**3
                    * const.h.cgs
                    / const.c.cgs**2
                )
                * (1 - nu_is / self.nu)
                * mean_intensities
            )

            integrated_heating_coefficient = integrate_array_by_blocks(
                heating_coefficient.values,
                self.nu,
                self.photoionization_block_references,
            )

        spontaneous_recombination_cooling_coefficient = (
            (
                8
                * np.pi
                * self.photoionization_cross_sections["x_sect"]
                * self.nu**3
                * const.h.cgs
                / const.c.cgs**2
            )
            * (1 - nu_i / self.nu)
            * boltzmann_factor
        )

        integrated_cooling_coefficient = integrate_array_by_blocks(
            spontaneous_recombination_cooling_coefficient.values,
            self.nu,
            self.photoionization_block_references,
        )

        heating_rate = integrated_heating_coefficient * level_number_density

        spontaneous_recombination_cooling_rate = (
            integrated_cooling_coefficient
            * phi
            * electron_density
            * ion_number_density
        )

        if stimulated_recombination_estimator is not None:
            stimulated_recombination_cooling_rate = (
                stimulated_recombination_estimator
                * phi
                * electron_density
                * ion_number_density
            )
        else:
            stimulated_recombination_cooling_rate = np.zeros(1)

        cooling_rate = (
            spontaneous_recombination_cooling_rate.sum()
            + stimulated_recombination_cooling_rate.sum()
        )

        return heating_rate, cooling_rate


class FreeFreeThermalRates:
    def __init__(self):
        # TODO replace with real value
        self.cooling_constant = 0

    def solve(self, heating_estimator, electron_temperature, heating_factor):
        if heating_estimator is not None:
            heating_rate = (
                heating_estimator
                * np.sqrt(electron_temperature)
                * heating_factor
            )
        else:
            heating_rate = 0

        cooling_rate = (
            self.cooling_constant
            * np.sqrt(electron_temperature)
            * heating_factor
        )

        return heating_rate, cooling_rate


class CollisionalIonizationThermalRates:
    def __init__(self, photoionization_cross_sections):
        self.nu = photoionization_cross_sections["nu"].values

    def solve(
        self,
        electron_density,
        collisional_ionization_rate_coefficient,
        phi,
        ion_number_density,
        level_number_density,
    ):
        heating_rate = (
            electron_density**2
            * collisional_ionization_rate_coefficient
            * phi
            * ion_number_density
            * self.nu
            * const.h.cgs
        )

        cooling_rate = (
            level_number_density
            * electron_density
            * collisional_ionization_rate_coefficient
        )

        return heating_rate, cooling_rate


class CollisionalBoundThermalRates:
    def __init__(self, collisional_cross_sections):
        self.nu = collisional_cross_sections["nu"].values

    def solve(
        self,
        electron_density,
        collisional_deexcitation_rate_coefficient,
        collisional_excitation_rate_coefficient,
        level_number_density,
    ):
        lower_level_number_density = level_number_density
        upper_level_number_density = level_number_density  # do something here

        heating_rate = (
            electron_density
            * collisional_deexcitation_rate_coefficient
            * upper_level_number_density
            * self.nu
            * const.h.cgs
        )

        cooling_rate = (
            collisional_excitation_rate_coefficient
            * self.nu
            * const.h.cgs
            * electron_density
            * lower_level_number_density
        )

        return heating_rate, cooling_rate


class AdiabaticThermalRates:
    """
    Class to represent the adiabatic heating and cooling rates.
    """

    def __init__(self):
        pass

    def solve(self, thermal_electron_distribution, time):
        return (
            3
            * const.k_B.cgs
            * thermal_electron_distribution.number_density
            * thermal_electron_distribution.temperature
            / time
        )
