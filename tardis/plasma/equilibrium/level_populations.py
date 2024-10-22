import numpy as np
import pandas as pd


class LevelPopulationSolver:
    def __init__(self, rates_matrices: pd.DataFrame, levels):
        self.rates_matrices = rates_matrices
        self.levels = levels

    def __calculate_boltzmann_factor(self, rates_matrix, species_id):
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
        normalized_level_boltzmann_factors = pd.DataFrame(
            index=self.levels.index, columns=self.rates_matrices.columns
        )

        for species_id in self.rates_matrices.index:
            # need to expand the numpy array so that it matches the dimensions of the normalized level boltzmann factors
            normalized_level_boltzmann_factors.loc[species_id, :] = (
                self.rates_matrices.loc[
                    species_id
                ].apply(
                    lambda rates_matrix: self.__calculate_boltzmann_factor(
                        rates_matrix, species_id
                    )
                )
            )

        return normalized_level_boltzmann_factors
