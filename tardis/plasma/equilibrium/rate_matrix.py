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
        thermal_electron_energy_distribution,
    ):
        """Construct the compiled rate matrix dataframe.

        Parameters
        ----------
        radiation_field : RadiationField
            Radiation field containing radiative temperature.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Distribution of electrons in the plasma, containing electron energies,
            temperatures and number densities.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell.
        """
        required_arg = {
            "radiative": radiation_field,
            "electron": thermal_electron_energy_distribution.temperature,
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

        # Multiply rates by electron number density where appropriate
        rates_df_list = [
            rates_df * thermal_electron_energy_distribution.number_density
            if solver_arg_tuple[1] == "electron"
            else rates_df
            for solver_arg_tuple, rates_df in zip(
                self.rate_solvers, rates_df_list
            )
        ]

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

        rate_matrices.index.names = ["atomic_number", "ion_number"]

        return rate_matrices


class IonRateMatrix:
    def __init__(
        self,
        radiative_ionization_rate_solver,
        collisional_ionization_rate_solver,
    ):
        """Constructs the ionization rate matrix from radiative and collisional
        ionization rate solvers.

        Parameters
        ----------
        radiative_ionization_rate_solver : AnalyticPhotoionizationRateSolver | EstimatedPhotoionizationRateSolver
            Solver for radiative ionization and recombination rates.
        collisional_ionization_rate_solver : CollisionalIonizationRateSolver
            Solver for collisional ionization and recombination rates.
        """
        self.radiative_ionization_rate_solver = radiative_ionization_rate_solver
        self.collisional_ionization_rate_solver = (
            collisional_ionization_rate_solver
        )

    def __calculate_total_grouped_rates(self, rates_df):
        """Helper function to calculate the total rates from the
        photoionization and recombination rates.

        Parameters
        ----------
        rates_df : pd.DataFrame
            DataFrame of rates indexed by atomic number and ion number,
            with each column being a cell.

        Returns
        -------
        pd.DataFrame
            A DataFrame of grouped total rates indexed by atomic number and ion number,
            with each column being a cell.
        """
        return (
            rates_df.groupby(
                level=(
                    "atomic_number",
                    "ion_number",
                    "ion_number_source",
                    "ion_number_destination",
                )
            )
            .sum()
            .groupby(level=("atomic_number"))
        )

    def __construct_rate_matrix(self, rate, cell, ion_states):
        """Construct a sparse rate matrix from the rates.

        Parameters
        ----------
        rate : pd.DataFrame
            Rate DataFrame indexed by atomic number and ion number
        shell : int
            Cell index
        ion_states : int
            Number of ion states for the atomic number

        Returns
        -------
        coo_matrix
            A sparse matrix representing the ionization rate for the given cell.
        """
        return coo_matrix(
            (
                rate[cell],
                (
                    rate.index.get_level_values("ion_number_source"),
                    rate.index.get_level_values("ion_number_destination"),
                ),
            ),
            shape=(ion_states, ion_states),
        )

    def solve(
        self,
        radiation_field,
        thermal_electron_energy_distribution,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        charge_conservation=False,
    ):
        """Compute the ionization rate matrix.

        Parameters
        ----------
        radiation_field : RadiationField
            A radiation field that can compute its mean intensity.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        lte_level_population : pd.DataFrame
            LTE level number density. Columns are cells.
        level_population : pd.DataFrame
            Estimated level number density. Columns are cells.
        lte_ion_population : pd.DataFrame
            LTE ion number density. Columns are cells.
        ion_population : pd.DataFrame
            Estimated ion number density. Columns are cells.
        charge_conservation : bool, optional
            Whether to include a charge conservation row in the rate matrix.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell. Each entry is a numpy array.
        """
        photoion_rates_df, recomb_rates_df = (
            self.radiative_ionization_rate_solver.solve(
                radiation_field,
                thermal_electron_energy_distribution,
                lte_level_population,
                level_population,
                lte_ion_population,
                ion_population,
            )
        )

        saha_factor = lte_level_population / (
            lte_ion_population.values
            * thermal_electron_energy_distribution.number_density.value
        )

        collisional_ionization_rates_df, collision_recombination_rates_df = (
            self.collisional_ionization_rate_solver.solve(
                thermal_electron_energy_distribution, saha_factor
            )
        )

        grouped_photoion_rates_df = self.__calculate_total_grouped_rates(
            photoion_rates_df
        )
        grouped_recomb_rates_df = self.__calculate_total_grouped_rates(
            recomb_rates_df
        )

        grouped_collisional_ionization_rates_df = (
            self.__calculate_total_grouped_rates(
                collisional_ionization_rates_df
            )
        )
        grouped_collisional_recombination_rates_df = (
            self.__calculate_total_grouped_rates(
                collision_recombination_rates_df
            )
        )

        rate_matrices = pd.DataFrame(
            index=grouped_photoion_rates_df.groups.keys(),
            columns=photoion_rates_df.columns,
        )

        for (atomic_number, photoion_rates), (
            atomic_number,
            recomb_rates,
        ), (atomic_number, coll_ion_rates), (
            atomic_number,
            recomb_ion_rates,
        ) in zip(
            grouped_photoion_rates_df,
            grouped_recomb_rates_df,
            grouped_collisional_ionization_rates_df,
            grouped_collisional_recombination_rates_df,
        ):
            ion_states = atomic_number + 1
            for shell in range(len(photoion_rates.columns)):
                photoion_matrix = self.__construct_rate_matrix(
                    photoion_rates, shell, ion_states
                )
                recomb_matrix = self.__construct_rate_matrix(
                    recomb_rates, shell, ion_states
                )
                coll_ion_matrix = self.__construct_rate_matrix(
                    coll_ion_rates, shell, ion_states
                )
                coll_recomb_matrix = self.__construct_rate_matrix(
                    recomb_ion_rates, shell, ion_states
                )

                matrix_array = (
                    photoion_matrix
                    + recomb_matrix
                    + coll_ion_matrix
                    + coll_recomb_matrix
                ).toarray()
                np.fill_diagonal(matrix_array, -np.sum(matrix_array, axis=0))
                matrix_array[0, :] = 1
                if charge_conservation:
                    charge_conservation_row = np.hstack(
                        (np.arange(0, ion_states), -1)
                    )
                    matrix_array = np.pad(matrix_array, ((0, 0), (0, 1)))
                    matrix_array = np.vstack(
                        (charge_conservation_row, matrix_array)
                    )
                rate_matrices.loc[atomic_number, shell] = matrix_array

        rate_matrices.index.names = ["atomic_number"]

        return rate_matrices
