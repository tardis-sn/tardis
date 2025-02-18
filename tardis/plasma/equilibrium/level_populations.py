import numpy as np
import pandas as pd


class LevelPopulationSolver:
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

    def __calculate_level_population(self, rates_matrix: np.ndarray):
        """Helper function to calculate the normalized, per-level boltzmann factor.

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
            rates_matrix[:, :], normalized_ion_population
        )
        return normalized_level_population

    def solve(self):
        """Solves the normalized level population values from the rate matrices.

        Returns
        -------
        pd.DataFrame
            Normalized level population values indexed by atomic number, ion
            number and level number. Columns are cells.
        """
        normalized_level_populations = pd.DataFrame(
            index=self.levels.index,
            columns=self.rates_matrices.columns,
            dtype=np.float64,
        )

        # try converting the set of vectors into a single 2D array and then applying index
        for species_id in self.rates_matrices.index:
            # TODO: resolve the chained assignment here. Maybe an intermediate df
            # is needed

            solved_matrices = self.rates_matrices.loc[species_id].apply(
                lambda rates_matrix: self.__calculate_level_population(
                    rates_matrix
                )
            )
            normalized_level_populations.loc[species_id, :] = np.vstack(
                solved_matrices.values
            ).T

        return normalized_level_populations
