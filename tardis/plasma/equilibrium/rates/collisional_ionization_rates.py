import pandas as pd

from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)
from tardis.plasma.equilibrium.rates.util import (
    align_ion_population_to_level_population,
    reindex_ionization_rate_dataframe,
)


class CollisionalIonizationRateSolver:
    """Solver for collisional ionization and recombination rates."""

    def __init__(self, photoionization_cross_sections):
        """Initialize the collisional ionization rate solver.

        Parameters
        ----------
        photoionization_cross_sections : pd.DataFrame
            Photoionization cross sections.
        """
        self.photoionization_cross_sections = photoionization_cross_sections

    def solve(
        self,
        electron_distribution,
        level_to_ion_population_factor,
        partition_function,
        level_boltzmann_factor,
        approximation="seaton",
    ):
        """Solve the collisional ionization and recombination rates.

        Parameters
        ----------
        electron_distribution : ThermalElectronEnergyDistribution
            Electron distribution per cell
        level_to_ion_population_factor : pandas.DataFrame, dtype float
            The level to ion population factor for each cell, Lucy 2003 Eq 14.
            Indexed by atom number, ion number, level number.
        approximation : str, optional
            The rate approximation to use, by default "seaton"

        Returns
        -------
        pd.DataFrame
            Collisional ionization rates
        pd.DataFrame
            Collisional recombination rates

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

        if isinstance(partition_function, pd.DataFrame):
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
