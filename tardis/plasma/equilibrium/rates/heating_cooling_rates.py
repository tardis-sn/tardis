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
        level_population: pd.DataFrame,
        ion_population: pd.DataFrame,
        thermal_electron_distribution: ThermalElectronEnergyDistribution,
        saha_factor: pd.DataFrame,
        radiation_field=None,
        bound_free_heating_estimator: pd.DataFrame | None = None,
        stimulated_recombination_estimator: pd.DataFrame | None = None,
    ) -> tuple[pd.Series, pd.Series]:
        """Compute the bound-free heating and cooling rates.

        Parameters
        ----------
        level_population : pd.DataFrame
            Estimated level number density. Columns represent cells.
        ion_population : pd.DataFrame
            Estimated ion number density. Columns represent cells.
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution containing the number density, temperature and energy.
        saha_factor : pd.DataFrame
            Saha factor for the ion populations as defined in Lucy 03 equation 14. Columns represent cells.
        radiation_field : RadiationField, optional
            A radiation field that can compute its mean intensity.
        bound_free_heating_estimator : pd.DataFrame, optional
            Montecarlo bound free heating estimator. Columns represent cells, by default None
        stimulated_recombination_estimator : pd.DataFrame, optional
            Montecarlo stimulated recombination estimator. Columns represent cells, by default None

        Returns
        -------
        tuple[pd.Series, pd.Series]
            Heating and cooling rates for the bound-free process for all cells.
        """
        nu_i = self.nu.groupby(level=[0, 1, 2]).first()
        nu_is = nu_i.loc[self.photoionization_cross_sections.index]
        n_cells = level_population.columns

        if bound_free_heating_estimator is not None:
            # TODO: check if this is correct
            integrated_heating_coefficient = bound_free_heating_estimator
        elif radiation_field is not None:
            mean_intensities = radiation_field.calculate_mean_intensity(
                self.nu.values * u.Hz
            )

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
            )
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

        integrated_cooling_coefficient = pd.DataFrame(
            integrate_array_by_blocks(
                spontaneous_recombination_cooling_coefficient.values,
                self.nu.values,
                self.photoionization_block_references,
            ),
            index=self.photoionization_index,
        )

        ion_cooling_factor = (
            thermal_electron_distribution.number_density.value
            * ion_population.loc[(1, 1)].values
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
            stimulated_recombination_cooling_rate = pd.DataFrame(
                np.zeros(
                    (len(spontaneous_recombination_cooling_rate), len(n_cells))
                ),
                index=spontaneous_recombination_cooling_rate.index,
                columns=n_cells,
            )

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
        self, ion_population: pd.DataFrame, electron_density: float
    ) -> pd.Series:
        """Compute the free-free heating factor.

        Parameters
        ----------
        ion_population : pd.DataFrame
            Ion number density. Columns represent cells.
        electron_density : float
            Electron number density value.

        Returns
        -------
        pd.Series
            The free-free heating factor for all cells.
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
        ion_population: pd.DataFrame,
    ) -> tuple[pd.Series, pd.Series]:
        """Compute the free-free heating and cooling rates for the input plasma conditions.

        Parameters
        ----------
        heating_estimator : float
            Montecarlo free-free heating estimator value.
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution containing the number density, temperature and energy.
        ion_population : pd.DataFrame
            Ion number density. Columns represent cells.

        Returns
        -------
        tuple[pd.Series, pd.Series]
            The heating and cooling rates for the free-free process for all cells.
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
        ion_population: pd.DataFrame,
        level_population: pd.DataFrame,
        collisional_ionization_rate_coefficient: pd.DataFrame,
        saha_factor: pd.DataFrame,
    ) -> tuple[pd.Series, pd.Series]:
        """Compute the collisional ionization heating and cooling rates.

        Parameters
        ----------
        electron_density : u.Quantity
            Electron number density with units.
        ion_population : pd.DataFrame
            Ion number density. Columns represent cells.
        level_population : pd.DataFrame
            Level number density. Columns represent cells.
        collisional_ionization_rate_coefficient : pd.DataFrame
            Collisional ionization rate coefficients. Columns represent cells.
        saha_factor : pd.DataFrame
            Saha factor for the ion populations as defined in Lucy 03 equation 14. Columns represent cells.

        Returns
        -------
        tuple[pd.Series, pd.Series]
            Heating and cooling rates for the collisional ionization process for all cells.
        """
        rate_factor = (
            electron_density.cgs.value
            * collisional_ionization_rate_coefficient.multiply(
                self.nu_i, axis=0
            )
            * const.h.cgs.value
        )

        heating_rate = (
            electron_density.cgs.value
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
        collisional_deexcitation_rate_coefficient: pd.DataFrame,
        collisional_excitation_rate_coefficient: pd.DataFrame,
        level_population: pd.DataFrame,
    ) -> tuple[pd.Series, pd.Series]:
        """Compute the collisional bound heating and cooling rates.

        Parameters
        ----------
        electron_density : u.Quantity
            Electron number density with units.
        collisional_deexcitation_rate_coefficient : pd.DataFrame
            Collisional deexcitation rate coefficients. Columns represent cells.
        collisional_excitation_rate_coefficient : pd.DataFrame
            Collisional excitation rate coefficients. Columns represent cells.
        level_population : pd.DataFrame
            Level number density. Columns represent cells.

        Returns
        -------
        tuple[pd.Series, pd.Series]
            Heating and cooling rates for the collisional bound process for all cells.
        """
        lower_index = collisional_excitation_rate_coefficient.index.droplevel(
            "level_number_upper"
        )
        upper_index = collisional_excitation_rate_coefficient.index.droplevel(
            "level_number_lower"
        )

        lower_level_number_density = level_population.loc[lower_index]
        upper_level_number_density = level_population.loc[upper_index]

        # For now, let's handle this similar to the original approach
        # to maintain compatibility, we'll loop through cells
        n_cells = level_population.columns
        heating_rates = []
        cooling_rates = []

        for cell in n_cells:
            heating_cell = (
                electron_density
                * collisional_deexcitation_rate_coefficient[cell]
                * upper_level_number_density[cell].values
                * self.nu
                * const.h.cgs
            ).sum()

            cooling_cell = (
                electron_density
                * collisional_excitation_rate_coefficient[cell]
                * lower_level_number_density[cell].values
                * self.nu
                * const.h.cgs
            ).sum()

            heating_rates.append(heating_cell)
            cooling_rates.append(cooling_cell)

        heating_rate = pd.Series(heating_rates, index=n_cells)
        cooling_rate = pd.Series(cooling_rates, index=n_cells)

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
