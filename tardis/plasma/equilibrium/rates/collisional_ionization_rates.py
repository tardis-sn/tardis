import astropy.units as u
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
        electron_distribution: ThermalElectronEnergyDistribution,
        level_to_ion_population_factor: pd.DataFrame,
        partition_function: pd.DataFrame | float,
        level_boltzmann_factor: pd.DataFrame,
        level_population: pd.DataFrame | None = None,
        ion_population: pd.DataFrame | None = None,
        approximation: str = "seaton",
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Solve the collisional ionization and recombination rates.

        Parameters
        ----------
        electron_distribution : ThermalElectronEnergyDistribution
            Electron distribution per cell
        level_to_ion_population_factor : pandas.DataFrame, dtype float
            The level to ion population factor for each cell, Lucy 2003 Eq 14.
            Indexed by atom number, ion number, level number.
        partition_function : pandas.DataFrame or float
            Partition function used for the LTE level fractions.
        level_boltzmann_factor : pandas.DataFrame
            Boltzmann factor used for the LTE level fractions.
        level_population : pandas.DataFrame, optional
            Estimated level populations used instead of LTE fractions.
        ion_population : pandas.DataFrame, optional
            Estimated ion populations used to normalize level populations.
        approximation : str, optional
            The rate approximation to use, by default "seaton"

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            Collisional ionization rates
            and collisional recombination rates.

        Raises
        ------
        ValueError
            If an unsupported approximation is requested.
        """
        collision_ionization_rates = self.solve_coefficients(
            electron_distribution.temperature, approximation
        )
        collision_ionization_rates.columns = (
            level_to_ion_population_factor.columns
        )

        # Inverse of the ionization rate for equilibrium
        collision_recombination_rates = collision_ionization_rates.multiply(
            level_to_ion_population_factor
        )

        if level_population is not None and ion_population is not None:
            level_population_fraction = level_population / (
                align_ion_population_to_level_population(
                    ion_population, level_population, next_higher=False
                )
            )
            level_population_fraction = level_population_fraction.loc[
                collision_ionization_rates.index
            ]
        else:
            if isinstance(partition_function, pd.DataFrame):
                partition_function = align_ion_population_to_level_population(
                    partition_function,
                    level_boltzmann_factor,
                    next_higher=False,
                )
            level_population_fraction = (
                level_boltzmann_factor / partition_function
            )

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

    def solve_coefficients(
        self, electron_temperature: u.Quantity, approximation: str = "seaton"
    ) -> pd.DataFrame:
        """Solve raw collisional-ionization rate coefficients."""
        if approximation != "seaton":
            raise ValueError(f"approximation {approximation} not supported")
        return CollisionalIonizationSeaton(
            self.photoionization_cross_sections
        ).solve(electron_temperature)
