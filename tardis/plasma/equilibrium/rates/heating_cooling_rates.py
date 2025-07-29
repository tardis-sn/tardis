import astropy.units as u
import numpy as np
import pandas as pd

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
        self.nu = photoionization_cross_sections["nu"]
        self.photoionization_block_references = np.pad(
            self.photoionization_cross_sections.nu.groupby(level=[0, 1, 2])
            .count()
            .values.cumsum(),
            [1, 0],
        )

    def solve(
        self,
        level_population,
        ion_population,
        thermal_electron_distribution,
        saha_factor,
        radiation_field=None,
        bound_free_heating_estimator=None,
        stimulated_recombination_estimator=None,
    ):
        """Compute the bound-free heating and cooling rates.

        Parameters
        ----------
        level_population : pd.DataFrame
            Estimated level number density. Columns are cells.
        ion_population : pd.DataFrame
            Estimated ion number density. Columns are cells.
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution containing the number density, temperature and energy.
        saha_factor : pd.DataFrame
            Saha factor for the ion populations as defined in Lucy 03 equation 14. Columns are cells.
        radiation_field : RadiationField, optional.
            A radiation field that can compute its mean intensity.
        bound_free_heating_estimator : pd.DataFrame, optional
            Montecarlo bound free heating estimator, by default None
        stimulated_recombination_estimator : pd.DataFrame, optional
            Montecarlo stimulated recombination estimator, by default None

        Returns
        -------
        pd.DataFrame, pd.DataFrame
            Heating and cooling rates for the bound-free process.
        """
        nu_i = (
            self.photoionization_cross_sections["nu"]
            .groupby(level=[0, 1, 2])
            .first()
        )
        nu_is = nu_i.loc[self.photoionization_cross_sections.index]

        if bound_free_heating_estimator is not None:
            # TODO: check if this is correct
            integrated_heating_coefficient = bound_free_heating_estimator
        else:
            mean_intensities = radiation_field.calculate_mean_intensity(self.nu)

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

            integrated_heating_coefficient = pd.DataFrame(
                integrate_array_by_blocks(
                    heating_coefficient.values,
                    self.nu.values,
                    self.photoionization_block_references,
                ),
                index=heating_coefficient.index,
                columns=heating_coefficient.columns,
            )

        boltzmann_factor = np.exp(
            -self.nu
            * u.Hz
            / thermal_electron_distribution.temperature
            * (const.h.cgs / const.k_B.cgs)
        )

        spontaneous_recombination_cooling_coefficient = pd.DataFrame(
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

        spontaneous_recombination_cooling_coefficient.insert(0, "nu", self.nu)

        integrated_cooling_coefficient = (
            spontaneous_recombination_cooling_coefficient.groupby(
                level=[0, 1, 2]
            ).apply(lambda sub: np.trapezoid(sub[0], sub["nu"]))
        )

        heating_rate = (
            integrated_heating_coefficient
            * level_population.loc[integrated_heating_coefficient.index]
        ).sum()

        spontaneous_recombination_cooling_rate = (
            integrated_cooling_coefficient
            * saha_factor.loc[integrated_cooling_coefficient.index]
            * thermal_electron_distribution.number_density.value
            * ion_population.loc[(1, 1)]  # Hydrogen ion population
        )

        if stimulated_recombination_estimator is not None:
            stimulated_recombination_cooling_rate = (
                stimulated_recombination_estimator
                * saha_factor.loc[stimulated_recombination_estimator.index]
                * thermal_electron_distribution.number_density.value
                * ion_population.loc[(1, 1)]
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
        self.cooling_constant = 1.426e-27  # in cgs units (see Osterbrock 1974)

    def heating_factor(self, ion_population, electron_density):
        """Compute the free-free heating factor.

        Parameters
        ----------
        ion_population : pd.DataFrame
            Ion number density. Columsn are cells.
        electron_density : pd.Series
            Electron number density.

        Returns
        -------
        pd.DataFrame
            The free-free heating factor.
        """
        ionic_charge_squared = np.square(
            ion_population.index.get_level_values(1).values
        )
        factor = (
            electron_density
            * ion_population.multiply(ionic_charge_squared, axis=0).sum()
        )
        return factor

    def solve(
        self,
        heating_estimator,
        thermal_electron_distribution,
        ion_population,
    ):
        """_summary_

        Parameters
        ----------
        heating_estimator : pd.DataFrame
            Montecarlo free-free heating estimator.
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution containing the number density, temperature and energy.
        ion_population : pd.DataFrame
            Ion number density. Columns are cells.

        Returns
        -------
        pd.DataFrame, pd.DataFrame
            The heating and cooling rates for the free-free process.
        """
        heating_factor = self.heating_factor(
            ion_population,
            thermal_electron_distribution.number_density,
        )

        heating_rate = (
            heating_estimator
            / np.sqrt(thermal_electron_distribution.temperature)
            * heating_factor
        )

        cooling_rate = (
            self.cooling_constant
            * np.sqrt(thermal_electron_distribution.temperature)
            * heating_factor
        )

        return heating_rate.value, cooling_rate.value


class CollisionalIonizationThermalRates:
    def __init__(self, photoionization_cross_sections):
        self.nu_i = (
            photoionization_cross_sections["nu"]
            .groupby(level=[0, 1, 2])
            .first()
        )

    def solve(
        self,
        electron_density,
        ion_population,
        level_population,
        collisional_ionization_rate_coefficient,
        saha_factor,
    ):
        """Compute the collisional ionization heating and cooling rates.

        Parameters
        ----------
        electron_density : pd.Series
            Electron number density.
        ion_population : pd.DataFrame
            Ion number density.
        level_population : pd.DataFrame
            Level number density.
        collisional_ionization_rate_coefficient : pd.DataFrame
            Collisional ionization rate ceofficients.
        saha_factor : pd.DataFrame
            Saha factor for the ion populations as defined in Lucy 03 equation 14. Columns are cells.

        Returns
        -------
        pd.DataFrame, pd.DataFrame
            Heating and cooling rates for the collisional ionization process.
        """
        rate_factor = (
            electron_density
            * collisional_ionization_rate_coefficient
            * self.nu_i
            * const.h.cgs.value
        )

        heating_rate = (
            electron_density
            * ion_population.loc[(1, 1)]
            * saha_factor
            * rate_factor
        ).sum()

        cooling_rate = (
            level_population.loc[collisional_ionization_rate_coefficient.index]
            * rate_factor
        ).sum()

        return heating_rate, cooling_rate


class CollisionalBoundThermalRates:
    def __init__(self, collisional_cross_sections):
        self.nu = collisional_cross_sections["nu"].values

    def solve(
        self,
        electron_density,
        collisional_deexcitation_rate_coefficient,
        collisional_excitation_rate_coefficient,
        level_population,
    ):
        """Compute the collisional bound heating and cooling rates.

        Parameters
        ----------
        electron_density : pd.Series
            Electron number density.
        collisional_deexcitation_rate_coefficient : pd.DataFrame
            Collisional deexcitation rate coefficients.
        collisional_excitation_rate_coefficient : pd.DataFrame
            Collisional excitation rate coefficients.
        level_population : pd.DataFrame
            Level number density. Columns are cells.

        Returns
        -------
        pd.DataFrame, pd.DataFrame
            Heating and cooling rates for the collisional bound process.
        """
        lower_index = collisional_excitation_rate_coefficient.index.droplevel(
            "level_number_upper"
        )
        upper_index = collisional_excitation_rate_coefficient.index.droplevel(
            "level_number_lower"
        )

        lower_level_number_density = level_population.loc[lower_index]
        upper_level_number_density = level_population.loc[upper_index]

        heating_rate = (
            electron_density
            * collisional_deexcitation_rate_coefficient
            * upper_level_number_density.values
            * self.nu
            * const.h.cgs
        ).sum()

        cooling_rate = (
            electron_density
            * collisional_excitation_rate_coefficient
            * lower_level_number_density.values
            * self.nu
            * const.h.cgs
        ).sum()

        return heating_rate, cooling_rate


class AdiabaticThermalRates:
    """
    Class to represent the adiabatic cooling rate.
    """

    def __init__(self):
        pass

    def solve(self, thermal_electron_distribution, time):
        """Solve for the adiabatic cooling rate.

        Parameters
        ----------
        thermal_electron_distribution : ThermalElectronDistribution
            The thermal electron distribution containing the number density and temperature.
        time : Quantity
            The time over which the adiabatic cooling is calculated.

        Returns
        -------
        Quantity
            The adiabatic cooling rate in TODO: add unit.
        """
        return (
            3
            * const.k_B.cgs
            * thermal_electron_distribution.number_density
            * thermal_electron_distribution.temperature
            / time
        )
