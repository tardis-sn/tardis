from __future__ import annotations

import astropy.units as u
import numpy as np
import pandas as pd

from tardis import constants as const
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.transport.montecarlo.estimators.util import (
    integrate_array_by_blocks,
)


class BoundFreeThermalRates:
    """
    Class to represent the bound-free heating rate.
    """

    def __init__(self, photoionization_cross_sections: pd.DataFrame) -> None:
        self.photoionization_cross_sections = photoionization_cross_sections
        self.nu = photoionization_cross_sections.nu
        self.photoionization_block_references = np.pad(
            self.photoionization_cross_sections.nu.groupby(level=[0, 1, 2])
            .count()
            .values.cumsum(),
            [1, 0],
        )

        self.photoionization_index = (
            self.photoionization_cross_sections.index.unique()
        )

    def solve(
        self,
        level_population: pd.Series,
        ion_population: pd.Series,
        thermal_electron_distribution: ThermalElectronEnergyDistribution,
        saha_factor: pd.Series,
        radiation_field=None,
        bound_free_heating_estimator: pd.Series | None = None,
        stimulated_recombination_estimator: pd.Series | None = None,
    ) -> tuple[float, float]:
        """Compute the bound-free heating and cooling rates.

        Parameters
        ----------
        level_population : pd.Series
            Estimated level number density for a single cell.
        ion_population : pd.Series
            Estimated ion number density for a single cell.
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution containing the number density, temperature and energy.
        saha_factor : pd.Series
            Saha factor for the ion populations as defined in Lucy 03 equation 14 for a single cell.
        radiation_field : RadiationField, optional
            A radiation field that can compute its mean intensity.
        bound_free_heating_estimator : pd.Series, optional
            Montecarlo bound free heating estimator for a single cell, by default None
        stimulated_recombination_estimator : pd.Series, optional
            Montecarlo stimulated recombination estimator for a single cell, by default None

        Returns
        -------
        tuple[float, float]
            Heating and cooling rates for the bound-free process.
        """
        nu_i = self.nu.groupby(level=[0, 1, 2]).first()
        nu_is = nu_i.loc[self.photoionization_cross_sections.index]

        if bound_free_heating_estimator is not None:
            # TODO: check if this is correct
            integrated_heating_coefficient = bound_free_heating_estimator
        elif radiation_field is not None:
            mean_intensities = radiation_field.calculate_mean_intensity(
                self.nu.values * u.Hz
            )[:, 0]

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
                    heating_coefficient.values[:, np.newaxis],
                    self.nu.values,
                    self.photoionization_block_references,
                ),
                index=self.photoionization_index,
            )[0]
        else:
            raise ValueError(
                "Either bound_free_heating_estimator or radiation_field must be provided."
            )

        boltzmann_factor = np.exp(
            -self.nu.values
            / thermal_electron_distribution.temperature.value
            * (const.h.cgs.value / const.k_B.cgs.value)
        )

        spontaneous_recombination_cooling_coefficient = pd.DataFrame(
            (
                8
                * np.pi
                * self.photoionization_cross_sections["x_sect"]
                * self.nu**3
                * const.h.cgs.value
                / const.c.cgs.value**2
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

        ion_cooling_factor = (
            thermal_electron_distribution.number_density.value
            * ion_population.loc[(1, 1)]
        )  # Hydrogen ion population

        spontaneous_recombination_cooling_rate = (
            integrated_cooling_coefficient
            * saha_factor.loc[integrated_cooling_coefficient.index]
            * ion_cooling_factor  # Hydrogen ion population
        )

        if stimulated_recombination_estimator is not None:
            stimulated_recombination_cooling_rate = (
                stimulated_recombination_estimator
                * saha_factor.loc[stimulated_recombination_estimator.index]
                * ion_cooling_factor
            )
        else:
            stimulated_recombination_cooling_rate = np.zeros(1)

        heating_rate = (
            integrated_heating_coefficient
            * level_population.loc[integrated_heating_coefficient.index]
        ).sum()

        cooling_rate = (
            spontaneous_recombination_cooling_rate
            + stimulated_recombination_cooling_rate
        ).sum()

        return heating_rate, cooling_rate


class FreeFreeThermalRates:
    def __init__(self) -> None:
        self.cooling_constant = 1.426e-27  # in cgs units (see Osterbrock 1974)

    def heating_factor(
        self, ion_population: pd.Series, electron_density: float
    ) -> pd.Series:
        """Compute the free-free heating factor.

        Parameters
        ----------
        ion_population : pd.Series
            Ion number density for a single cell.
        electron_density : float
            Electron number density value.

        Returns
        -------
        pd.Series
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
        heating_estimator: float,
        thermal_electron_distribution: ThermalElectronEnergyDistribution,
        ion_population: pd.Series,
    ) -> tuple[float, float]:
        """Compute the free-free heating and cooling rates for the input plasma conditions.

        Parameters
        ----------
        heating_estimator : float
            Montecarlo free-free heating estimator value.
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution containing the number density, temperature and energy.
        ion_population : pd.Series
            Ion number density for a single cell.

        Returns
        -------
        tuple[float, float]
            The heating and cooling rates for the free-free process.
        """
        heating_factor = self.heating_factor(
            ion_population,
            thermal_electron_distribution.number_density.cgs.value,
        )

        heating_rate = (
            heating_estimator
            / np.sqrt(thermal_electron_distribution.temperature.cgs.value)
            * heating_factor
        )

        cooling_rate = (
            self.cooling_constant
            * np.sqrt(thermal_electron_distribution.temperature.cgs.value)
            * heating_factor
        )

        return heating_rate, cooling_rate


class CollisionalIonizationThermalRates:
    def __init__(self, photoionization_cross_sections: pd.DataFrame) -> None:
        self.nu_i = (
            photoionization_cross_sections["nu"]
            .groupby(level=[0, 1, 2])
            .first()
        )

    def solve(
        self,
        electron_density: u.Quantity,
        ion_population: pd.Series,
        level_population: pd.Series,
        collisional_ionization_rate_coefficient: pd.Series,
        saha_factor: pd.Series,
    ) -> tuple[float, float]:
        """Compute the collisional ionization heating and cooling rates.

        Parameters
        ----------
        electron_density : u.Quantity
            Electron number density with units.
        ion_population : pd.Series
            Ion number density for a single cell.
        level_population : pd.Series
            Level number density for a single cell.
        collisional_ionization_rate_coefficient : pd.Series
            Collisional ionization rate coefficients for a single cell.
        saha_factor : pd.Series
            Saha factor for the ion populations as defined in Lucy 03 equation 14 for a single cell.

        Returns
        -------
        tuple[float, float]
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
    def __init__(self, lines: pd.DataFrame) -> None:
        self.nu = lines["nu"].values

    def solve(
        self,
        electron_density: u.Quantity,
        collisional_deexcitation_rate_coefficient: pd.Series,
        collisional_excitation_rate_coefficient: pd.Series,
        level_population: pd.Series,
    ) -> tuple[float, float]:
        """Compute the collisional bound heating and cooling rates.

        Parameters
        ----------
        electron_density : u.Quantity
            Electron number density with units.
        collisional_deexcitation_rate_coefficient : pd.Series
            Collisional deexcitation rate coefficients for a single cell.
        collisional_excitation_rate_coefficient : pd.Series
            Collisional excitation rate coefficients for a single cell.
        level_population : pd.Series
            Level number density for a single cell.

        Returns
        -------
        tuple[float, float]
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

    def __init__(self) -> None:
        pass

    def solve(
        self,
        thermal_electron_distribution: ThermalElectronEnergyDistribution,
        time: u.Quantity,
    ) -> u.Quantity:
        """Solve for the adiabatic cooling rate.

        Parameters
        ----------
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            The thermal electron distribution containing the number density and temperature.
        time : u.Quantity
            The time over which the adiabatic cooling is calculated.

        Returns
        -------
        u.Quantity
            The adiabatic cooling rate in erg cm^-3 s^-1.
        """
        return (
            3
            * const.k_B.cgs
            * thermal_electron_distribution.number_density
            * thermal_electron_distribution.temperature
            / time
        )
