from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

    from tardis.plasma.electron_energy_distribution import (
        ThermalElectronEnergyDistribution,
    )
    from tardis.plasma.radiation_field import DilutePlanckianRadiationField


class ThermalBalanceSolver:
    """Class to solve the thermal balance equation using all available
    heating and cooling processes. See section 6.4 in Lucy 03.
    """

    def __init__(self, thermal_rates: dict) -> None:
        """Initialize with thermal rates.

        Parameters
        ----------
        thermal_rates : dict
            Thermal rate objects with keys: "bound_free", "free_free",
            "collisional_ionization", and "collisional_bound".
        """
        self.thermal_rates = thermal_rates

    def solve(
        self,
        thermal_electron_distribution: ThermalElectronEnergyDistribution,
        level_population: pd.DataFrame,
        ion_population: pd.DataFrame,
        collisional_ionization_rate_coefficient: pd.DataFrame,
        collisional_deexcitation_rate_coefficient: pd.DataFrame,
        collisional_excitation_rate_coefficient: pd.DataFrame,
        free_free_heating_estimator: pd.DataFrame,
        level_population_ratio: pd.DataFrame,
        radiation_field: DilutePlanckianRadiationField | None,
        bound_free_heating_estimator: pd.DataFrame | None = None,
        stimulated_recombination_estimator: pd.DataFrame | None = None,
    ) -> tuple[pd.Series, pd.Series]:
        """Compute the current heating rate and the fractional heating
        rate using all available processes. See section 6.4 in Lucy 03.

        Parameters
        ----------
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy, temperature, and density.
        level_population : pd.DataFrame
            Level number density.
        ion_population : pd.DataFrame
            Ion number density.
        collisional_ionization_rate_coefficient : pd.DataFrame
            Collisional ionization rate coefficient.
        collisional_deexcitation_rate_coefficient : pd.DataFrame
            Collisional deexcitation rate coefficient.
        collisional_excitation_rate_coefficient : pd.DataFrame
            Collisional excitation rate coefficient.
        free_free_heating_estimator : pd.DataFrame
            Montecarlo estimator for free-free heating.
        level_population_ratio : pd.DataFrame
            Level population to ion population ratio. Lucy 03, equation 14.
        radiation_field : RadiationField, optional
            Radiation field for mean intensity calculation.
        bound_free_heating_estimator : pd.DataFrame, optional
            Bound-free heating estimator.
        stimulated_recombination_estimator : pd.DataFrame, optional
            Stimulated recombination estimator.

        Returns
        -------
        tuple[pd.Series, pd.Series]
            Total heating rate and fractional heating rate for each cell.
        """
        electron_density = thermal_electron_distribution.number_density

        bound_free_heating, free_bound_cooling = self.thermal_rates[
            "bound_free"
        ].solve(
            level_population,
            ion_population,
            thermal_electron_distribution,
            level_population_ratio,
            radiation_field,
            bound_free_heating_estimator,
            stimulated_recombination_estimator,
        )

        free_free_heating, free_free_cooling = self.thermal_rates[
            "free_free"
        ].solve(
            free_free_heating_estimator,
            thermal_electron_distribution,
            ion_population,
        )

        collisional_ionization_heating, collisional_ionization_cooling = (
            self.thermal_rates["collisional_ionization"].solve(
                electron_density,
                ion_population,
                level_population,
                collisional_ionization_rate_coefficient,
                level_population_ratio,
            )
        )

        collisional_bound_heating, collisional_bound_cooling = (
            self.thermal_rates["collisional_bound"].solve(
                electron_density,
                collisional_deexcitation_rate_coefficient,
                collisional_excitation_rate_coefficient,
                level_population,
            )
        )

        total_heating = (
            bound_free_heating
            + free_free_heating
            + collisional_ionization_heating
            + collisional_bound_heating
        )
        total_cooling = (
            free_bound_cooling
            + free_free_cooling
            + collisional_ionization_cooling
            + collisional_bound_cooling
        )

        total_heating_rate = total_heating - total_cooling
        fractional_heating_rate = (
            total_heating - total_cooling
        ) / total_cooling

        return total_heating_rate, fractional_heating_rate
