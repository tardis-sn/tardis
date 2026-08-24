import pandas as pd

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)
from tardis.plasma.equilibrium.rates.util import (
    align_ion_population_to_level_population,
    reindex_ionization_rate_dataframe,
)


class CollisionalIonizationRateSolver:
    """Solver for collisional ionization and recombination rates."""

    def __init__(
        self, photoionization_cross_sections: pd.DataFrame
    ) -> None:
        """Initialize the collisional ionization rate solver.

        Parameters
        ----------
        photoionization_cross_sections : pd.DataFrame
            Photoionization cross sections indexed by atomic number, ion
            number, and level number.
        """
        self.photoionization_cross_sections = photoionization_cross_sections

    def solve(
        self,
        electron_distribution: ThermalElectronEnergyDistribution,
        level_to_ion_population_factor: pd.DataFrame,
        partition_function: pd.DataFrame,
        level_boltzmann_factor: pd.DataFrame,
        approximation: str = "seaton",
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Solve the collisional ionization and recombination rates.

        Parameters
        ----------
        electron_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution per cell.
        level_to_ion_population_factor : pd.DataFrame
            The level to ion population factor for each cell, Lucy 2003 Eq 14.
            Indexed by atom number, ion number, level number.
        partition_function : pd.DataFrame
            Partition function for each ion and cell.
        level_boltzmann_factor : pd.DataFrame
            Boltzmann factor for each level and cell.
        approximation : str, optional
            The rate approximation to use, by default ``"seaton"``.

        Returns
        -------
        tuple[pd.DataFrame, pd.DataFrame]
            Collisional ionization rates and collisional recombination rates.

        Raises
        ------
        ValueError
            If an unsupported approximation is requested.
        """
        if approximation == "seaton":
            strength_solver = CollisionalIonizationSeaton(
                self.photoionization_cross_sections
            )
        else:
            raise ValueError(f"approximation {approximation} not supported")

        collision_ionization_rates = strength_solver.solve(
            electron_distribution.temperature
        )
        collision_ionization_rates.columns = (
            level_to_ion_population_factor.columns
        )

        # Inverse of the ionization rate for equilibrium
        collision_recombination_rates = collision_ionization_rates.multiply(
            level_to_ion_population_factor
        )

        partition_function = align_ion_population_to_level_population(
            partition_function,
            level_boltzmann_factor,
            next_higher=False,
        )
        level_population_fraction = level_boltzmann_factor / partition_function

        # used to scale the photoionization rate because we keep the level population
        # fixed while we calculated the ion number density
        collision_ionization_rates = (
            reindex_ionization_rate_dataframe(
                collision_ionization_rates * level_population_fraction,
                recombination=False,
            )
            * electron_distribution.number_density
        )

        collision_recombination_rates = (
            reindex_ionization_rate_dataframe(
                collision_recombination_rates, recombination=True
            )
        ) * electron_distribution.number_density**2

        return collision_ionization_rates, collision_recombination_rates
