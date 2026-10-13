import astropy.units as u
import pandas as pd

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_ion_number_density_to_level_number_density,
    reindex_ionization_rate_dataframe,
)


class CollisionalIonizationRateSolver:
    """Solver for collisional ionization and recombination rates."""

    def __init__(self, photoionization_cross_sections: pd.DataFrame) -> None:
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
        level_to_ion_number_density_factor: pd.DataFrame,
        partition_function: pd.DataFrame,
        level_boltzmann_factor: pd.DataFrame,
        level_number_density: pd.DataFrame | None = None,
        ion_number_density: pd.DataFrame | None = None,
        approximation: str = "seaton",
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Solve the collisional ionization and recombination rates.

        Parameters
        ----------
        electron_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution per cell.
        level_to_ion_number_density_factor : pd.DataFrame
            The level to ion number density factor for each cell, Lucy 2003 Eq 14.
            Indexed by atom number, ion number, level number.
        partition_function : pd.DataFrame
            Partition function for each ion and cell.
        level_boltzmann_factor : pd.DataFrame
            Boltzmann factor for each level and cell.
        level_number_density : pandas.DataFrame, optional
            Estimated level number densities used instead of LTE fractions.
        ion_number_density : pandas.DataFrame, optional
            Estimated ion number densities used to normalize level number densities.
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
        collision_ionization_rates = self.solve_coefficients(
            electron_distribution.temperature, approximation
        )
        collision_ionization_rates.columns = (
            level_to_ion_number_density_factor.columns
        )

        # Inverse of the ionization rate for equilibrium
        collision_recombination_rates = collision_ionization_rates.multiply(
            level_to_ion_number_density_factor
        )

        if level_number_density is not None and ion_number_density is not None:
            fractional_level_number_density = level_number_density / (
                reindex_ion_number_density_to_level_number_density(
                    ion_number_density, level_number_density, next_higher=False
                )
            )
            fractional_level_number_density = fractional_level_number_density.loc[
                collision_ionization_rates.index
            ]
        else:
            partition_function = reindex_ion_number_density_to_level_number_density(
                partition_function,
                level_boltzmann_factor,
                next_higher=False,
            )
            fractional_level_number_density = (
                level_boltzmann_factor / partition_function
            )

        # used to scale the photoionization rate because we keep the level number density
        # fixed while we calculated the ion number density
        collision_ionization_rates = (
            reindex_ionization_rate_dataframe(
                collision_ionization_rates * fractional_level_number_density,
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
