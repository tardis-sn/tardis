import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix


class RateMatrix:
    def __init__(
        self,
        rate_solvers: list,
        levels: pd.DataFrame,
    ):
        """Constructs the rate matrix from an arbitrary number of rate solvers.

        Parameters
        ----------
        rate_solvers : list
            List of rate solver objects.
        levels : pd.DataFrame
            DataFrame of energy levels.
        """
        self.rate_solvers = rate_solvers
        self.levels = levels

    def solve(
        self,
        radiation_field,
        electron_distribution,
    ):
        """Construct the compiled rate matrix dataframe.

        Parameters
        ----------
        radiation_field : RadiationField
            Radiation field containing radiative temperature.
        electron_distribution : ElectronDistribution
            Distribution of electrons in the plasma, containing electron
            temperatures and number densities.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell.
        """
        required_arg = {
            "radiative": radiation_field,
            "electron": electron_distribution.temperature,
        }

        rates_df_list = [
            solver.solve(required_arg[arg]) for solver, arg in self.rate_solvers
        ]
        # Extract all indexes
        all_indexes = set()
        for df in rates_df_list:
            all_indexes.update(df.index)

        # Create a union of all indexes
        all_indexes = sorted(all_indexes)

        # Reindex each dataframe to ensure consistent indices
        rates_df_list = [
            df.reindex(all_indexes, fill_value=0) for df in rates_df_list
        ]

        rates_df_list[1] *= electron_distribution.number_density
        rates_df = sum(rates_df_list)

        grouped_rates_df = rates_df.groupby(
            level=("atomic_number", "ion_number")
        )

        rate_matrices = pd.DataFrame(
            index=grouped_rates_df.groups.keys(), columns=rates_df.columns
        )

        for species_id, rates in grouped_rates_df:
            number_of_levels = self.levels.energy.loc[species_id].count()
            for shell in range(len(rates.columns)):
                matrix = coo_matrix(
                    (
                        rates[shell],
                        (
                            rates.index.get_level_values(
                                "level_number_destination"
                            ),
                            rates.index.get_level_values("level_number_source"),
                        ),
                    ),
                    shape=(number_of_levels, number_of_levels),
                )
                matrix_array = matrix.toarray()
                np.fill_diagonal(matrix_array, -np.sum(matrix_array, axis=0))
                matrix_array[0, :] = 1
                rate_matrices.loc[species_id, shell] = matrix_array

        return rate_matrices
