import numpy as np
import numpy.typing as npt
import pandas as pd


class LevelNumberDensitySolver:
    """Solve fractional level number densities from bound-bound rate matrices."""

    def __init__(self, rates_matrices: pd.DataFrame, levels: pd.DataFrame):
        """Solve the fractional level number density values from the rate matrices.

        Parameters
        ----------
        rates_matrices : pd.DataFrame
            DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell.
        levels : pd.DataFrame
            DataFrame of energy levels.
        """
        self.rates_matrices = rates_matrices
        self.levels = levels

    def _calculate_level_number_density(
        self, rates_matrix: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Calculate fractional per-level number densities.

        Parameters
        ----------
        rates_matrix : np.ndarray
            The rate matrix for a given species and cell.

        Returns
        -------
        np.ndarray
            The normalized, per-level number density.
        """
        fractional_ion_number_density = np.zeros(rates_matrix.shape[0])
        fractional_ion_number_density[0] = 1.0
        fractional_level_number_density = np.linalg.solve(
            rates_matrix, fractional_ion_number_density
        )
        return fractional_level_number_density

    def solve(self) -> pd.DataFrame:
        """Solves the fractional level number density values from the rate matrices.

        Returns
        -------
        pd.DataFrame
            Normalized level number density values indexed by atomic number, ion
            number and level number. Columns are cells.
        """
        fractional_level_number_densities = np.full(
            (len(self.levels), len(self.rates_matrices.columns)), np.nan
        )

        for species_id in self.rates_matrices.index:
            matrices = np.stack(
                self.rates_matrices.loc[species_id].to_numpy()
            )
            number_densities = np.array(
                [self._calculate_level_number_density(matrix) for matrix in matrices]
            ).T
            fractional_level_number_densities[
                self.levels.index.get_loc(species_id)
            ] = number_densities

        return pd.DataFrame(
            fractional_level_number_densities,
            index=self.levels.index,
            columns=self.rates_matrices.columns,
        )
