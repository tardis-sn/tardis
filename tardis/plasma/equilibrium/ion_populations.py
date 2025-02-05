import numpy as np
import pandas as pd


class IonPopulationSolver:
    def __init__(self, rates_matrices: pd.DataFrame, ions: pd.DataFrame):
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rates_matrices : pd.DataFrame
            DataFrame of rate matrices indexed by atomic number,
            with each column being a cell.
        ions : pd.DataFrame
            DataFrame of ions present.
        """
        self.rates_matrices = rates_matrices
        self.ions = ions

    def __calculate_ion_population(self, rates_matrix: np.ndarray):
        """Helper function to calculate the normalized, per-ion boltzmann factor.

        Parameters
        ----------
        rates_matrix : np.ndarray
            The rate matrix for a given species and cell.

        Returns
        -------
        np.ndarray
            The normalized, per-ion population.
        """
        normalized_ion_population = np.zeros(rates_matrix.shape[0])
        normalized_ion_population[0] = 1.0
        normalized_ion_population = np.linalg.solve(
            rates_matrix, normalized_ion_population
        )
        return normalized_ion_population

    def solve(self):
        """Solves the normalized ion population values from the rate matrices.

        Returns
        -------
        pd.DataFrame
            Normalized ion population values indexed by atomic number, ion
            number and ion number. Columns are cells.
        """
        normalized_ion_populations = pd.DataFrame(
            index=pd.MultiIndex.from_tuples(self.ions),
            columns=self.rates_matrices.columns,
            dtype=np.float64,
        )

        # try converting the set of vectors into a single 2D array and then applying index
        for species_id in self.rates_matrices.index:
            # TODO: resolve the chained assignment here. Maybe an intermediate df
            # is needed

            solved_matrices = self.rates_matrices.loc[species_id].apply(
                lambda rates_matrix: self.__calculate_ion_population(
                    rates_matrix
                )
            )
            normalized_ion_populations.loc[species_id, :] = np.vstack(
                solved_matrices.values
            ).T

        return normalized_ion_populations
