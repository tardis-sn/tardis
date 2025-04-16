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
    ):
        self.radiative_ionization_rate_solver = radiative_ionization_rate_solver

    def solve(
        self,
        radiation_field,
        thermal_electron_energy_distribution,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
    ):
        """Compute the ionization rate matrix in the
        case where the radiation field is not estimated.

        Parameters
        ----------
        radiation_field : RadiationField
            A radiation field that can compute its mean intensity.
        electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        level_number_density : pd.DataFrame
            Electron energy level number density. Columns are cells.
        ion_number_density : pd.DataFrame
            Ion number density. Columns are cells.
        saha_factor : pd.DataFrame
            Saha factor: the LTE level number density divided by the LTE ion
            number density and the electron number density.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell.
        """
        ion_rates_df, recomb_rates_df = (
            self.radiative_ionization_rate_solver.solve(
                radiation_field,
                thermal_electron_energy_distribution,
                lte_level_population,
                level_population,
                lte_ion_population,
                ion_population,
            )
        )

        ion_rates_total = ion_rates_df.groupby(
            level=(
                "atomic_number",
                "ion_number",
                "ion_number_source",
                "ion_number_destination",
            )
        ).sum()

        recomb_rates_total = recomb_rates_df.groupby(
            level=(
                "atomic_number",
                "ion_number",
                "ion_number_source",
                "ion_number_destination",
            )
        ).sum()

        grouped_ion_rates_df = ion_rates_total.groupby(level=("atomic_number"))
        grouped_recomb_rates_df = recomb_rates_total.groupby(
            level=("atomic_number")
        )

        rate_matrices = pd.DataFrame(
            index=grouped_ion_rates_df.groups.keys(),
            columns=ion_rates_df.columns,
        )

        for (atomic_number, ion_rates), (atomic_number, recomb_rates) in zip(
            grouped_ion_rates_df, grouped_recomb_rates_df
        ):
            ion_states = atomic_number + 1
            for shell in range(len(ion_rates.columns)):
                ion_matrix = coo_matrix(
                    (
                        ion_rates[shell],
                        (
                            ion_rates.index.get_level_values(
                                "ion_number_source"
                            ),
                            ion_rates.index.get_level_values(
                                "ion_number_destination"
                            ),
                        ),
                    ),
                    shape=(ion_states, ion_states + 1),
                )
                recomb_matrix = coo_matrix(
                    (
                        recomb_rates[shell],
                        (
                            recomb_rates.index.get_level_values(
                                "ion_number_source"
                            ),
                            recomb_rates.index.get_level_values(
                                "ion_number_destination"
                            ),
                        ),
                    ),
                    shape=(ion_states, ion_states + 1),
                )

                matrix_array = (ion_matrix + recomb_matrix).toarray()
                np.fill_diagonal(matrix_array, -np.sum(matrix_array, axis=0))
                charge_conservation_row = np.hstack(
                    (np.arange(0, ion_states), -1)
                )
                matrix_array = np.vstack(
                    (charge_conservation_row, matrix_array)
                )
                rate_matrices.loc[atomic_number, shell] = matrix_array

        rate_matrices.index.names = ["atomic_number"]

        return rate_matrices


class IonFullRateMatrix:
    def __init__(
        self,
        levels_rate_solvers,
        levels,
        ion_rate_solver,
        photoionization_cross_sections,
    ):
        self.levels_rate_solvers = levels_rate_solvers
        self.levels = levels
        self.ion_rate_solver = ion_rate_solver
        self.photoionization_cross_sections = photoionization_cross_sections

    def solve(
        self,
        radiation_field,
        thermal_electron_energy_distribution,
        level_number_density,
        ion_number_density,
        saha_factor,
    ):
        """Compute the ionization rate matrix in the
        case where the radiation field is not estimated.

        Parameters
        ----------
        radiation_field : RadiationField
            A radiation field that can compute its mean intensity.
        electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        level_number_density : pd.DataFrame
            Electron energy level number density. Columns are cells.
        ion_number_density : pd.DataFrame
            Ion number density. Columns are cells.
        saha_factor : pd.DataFrame
            Saha factor: the LTE level number density divided by the LTE ion
            number density and the electron number density.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell.
        """
        ion_rates_df = self.ion_rate_solver.solve(
            radiation_field,
            thermal_electron_energy_distribution,
            level_number_density,
            ion_number_density,
            saha_factor,
        )

        # Ionization rates R_ik connect level i to ion k, level i + j where j is
        # the number of levels in ion k-1. Thus photoionization rates can be summed
        # with radiative rates.
        # But we don't know which level in the next ion the electron goes to?
        # Seems that they could be treated as the LTE level number density, see
        # Mihalas 2014 eq 9.84 and discussion
        # Chris Fontes says not to bother if the data doesn't include dest. levels

        grouped_rates_df = ion_rates_df.groupby(
            level=("atomic_number", "ion_number")
        )

        rate_matrices = pd.DataFrame(
            index=grouped_rates_df.groups.keys(), columns=ion_rates_df.columns
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
