import numpy as np
import pandas as pd


def expand_series(series: pd.Series):
    """Helper function to expand a series of numpy arrays into a DataFrame.

    Parameters
    ----------
    series : pd.Series
        A series of numpy arrays.

    Returns
    -------
    pd.DataFrame
        The expanded DataFrame where each row is an entry from the numpy array.
    """
    # Expand the numpy arrays into rows
    expanded_series = pd.DataFrame(
        np.concatenate(series.values).reshape(-1, len(series))
    )
    return expanded_series


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

    def __calculate_boltzmann_factor(
        self, rates_matrix: np.ndarray, species_id: tuple
    ):
        """Helper function to calculate the normalized, per-level boltzmann factor.

        Parameters
        ----------
        rates_matrix : np.ndarray
            The rate matrix for a given species and cell.
        species_id : tuple
            A tuple of the atomic number and ion number.

        Returns
        -------
        np.ndarray
            The normalized, per-level boltzmann factor.
        """
        x = np.zeros(rates_matrix.shape[0])
        x[0] = 1.0
        level_boltzmann_factor = np.linalg.solve(rates_matrix[:, :], x)
        normalized_level_boltzmann_factor = (
            level_boltzmann_factor
            * self.levels.g.loc[species_id][0]
            / level_boltzmann_factor[0]
        )
        return normalized_level_boltzmann_factor

    def solve(self):
        """Solves the normalized level population values from the rate matrices.

        Returns
        -------
        pd.DataFrame
            Normalized level population values indexed by atomic number, ion
            number and level number. Columns are cells.
        """
        normalized_level_boltzmann_factors = pd.DataFrame(
            index=self.levels.index,
            columns=self.rates_matrices.columns,
            dtype=np.float64,
        )

        for species_id in self.rates_matrices.index:
            # TODO: resolve the chained assignment here. Maybe an intermediate df
            # is needed
            normalized_level_boltzmann_factors.loc[species_id, :].update(
                expand_series(
                    self.rates_matrices.loc[species_id].apply(
                        lambda rates_matrix: self.__calculate_boltzmann_factor(
                            rates_matrix, species_id
                        )
                    )
                )
            )

        return normalized_level_boltzmann_factors
