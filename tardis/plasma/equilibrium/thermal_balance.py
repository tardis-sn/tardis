import pandas as pd

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField


class ThermalBalanceSolver:
    """Class to solve the thermal balance equation using all available
    heating and cooling processes. See section 6.4 in Lucy 03, especially
    equations 55-57.
    """

    def __init__(
        self,
        bound_free_solver: BoundFreeThermalRates,
        free_free_solver: FreeFreeThermalRates,
        collisional_ionization_solver: CollisionalIonizationThermalRates,
        collisional_bound_solver: CollisionalBoundThermalRates,
    ) -> None:
        """Initialize the thermal balance solver with individual solver components.

        Parameters
        ----------
        bound_free_solver
            Solver for bound-free heating and free-bound cooling processes.
        free_free_solver
            Solver for free-free heating and cooling processes.
        collisional_ionization_solver
            Solver for collisional ionization heating and cooling processes.
        collisional_bound_solver
            Solver for collisional bound-bound heating and cooling processes.
        """
        self.bound_free_solver = bound_free_solver
        self.free_free_solver = free_free_solver
        self.collisional_ionization_solver = collisional_ionization_solver
        self.collisional_bound_solver = collisional_bound_solver

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
        radiation_field: DilutePlanckianRadiationField | None = None,
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

        bound_free_heating, free_bound_cooling = self.bound_free_solver.solve(
            level_population,
            ion_population,
            thermal_electron_distribution,
            level_population_ratio,
            radiation_field,
            bound_free_heating_estimator,
            stimulated_recombination_estimator,
        )

        free_free_heating, free_free_cooling = self.free_free_solver.solve(
            free_free_heating_estimator,
            thermal_electron_distribution,
            ion_population,
        )

        collisional_ionization_heating, collisional_ionization_cooling = (
            self.collisional_ionization_solver.solve(
                electron_density,
                ion_population,
                level_population,
                collisional_ionization_rate_coefficient,
                level_population_ratio,
            )
        )

        collisional_bound_heating, collisional_bound_cooling = (
            self.collisional_bound_solver.solve(
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
