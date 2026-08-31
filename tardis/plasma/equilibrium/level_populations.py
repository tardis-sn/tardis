import numpy as np
import numpy.typing as npt
import pandas as pd


class LevelPopulationSolver:
    """Solve normalized level populations from bound-bound rate matrices."""

    def __init__(self, rates_matrices: pd.DataFrame, levels: pd.DataFrame):
        """Solve the normalized level population values from the rate matrices.

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

    def _calculate_level_population(
        self, rates_matrix: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Calculate normalized per-level populations.

        Parameters
        ----------
        rates_matrix : np.ndarray
            The rate matrix for a given species and cell.

        Returns
        -------
        np.ndarray
            The normalized, per-level population.
        """
        normalized_ion_population = np.zeros(rates_matrix.shape[0])
        normalized_ion_population[0] = 1.0
        normalized_level_population = np.linalg.solve(
            rates_matrix, normalized_ion_population
        )
        return normalized_level_population

    def solve(self) -> pd.DataFrame:
        """Solves the normalized level population values from the rate matrices.

        Returns
        -------
        pd.DataFrame
            Normalized level population values indexed by atomic number, ion
            number and level number. Columns are cells.
        """
        normalized_level_populations = np.full(
            (len(self.levels), len(self.rates_matrices.columns)), np.nan
        )

        for species_id in self.rates_matrices.index:
            matrices = np.stack(
                self.rates_matrices.loc[species_id].to_numpy()
            )
            populations = np.array(
                [self._calculate_level_population(matrix) for matrix in matrices]
            ).T
            normalized_level_populations[
                self.levels.index.get_loc(species_id)
            ] = populations

        return pd.DataFrame(
            normalized_level_populations,
            index=self.levels.index,
            columns=self.rates_matrices.columns,
        )
