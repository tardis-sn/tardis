import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix


class LevelRateMatrix:
    def __init__(
        self,
        levels: pd.DataFrame,
    ):
        """Constructs the level rate matrix from precomputed rate data.

        Parameters
        ----------
        levels : pd.DataFrame
            DataFrame of energy levels.
        """
        self.levels = levels

    def solve(
        self,
        radiative_rates_df: pd.DataFrame,
        collisional_rates_df: pd.DataFrame,
        electron_number_density,
    ):
        """Construct the compiled rate matrix dataframe.

        Parameters
        ----------
        radiative_rates_df : pd.DataFrame
            Precomputed radiative rates indexed by transition identifiers,
            with each column being a cell.
        collisional_rates_df : pd.DataFrame
            Precomputed thermal collisional rates indexed by transition
            identifiers, with each column being a cell.
        electron_number_density : pandas.Series | numpy.ndarray
            Electron number density indexed by cell.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell.
        """
        rates_df_list = [radiative_rates_df, collisional_rates_df]
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
        rates_df_list[1] *= electron_number_density

        rates_df = rates_df_list[0].copy()
        for rates_df_component in rates_df_list[1:]:
            rates_df = rates_df.add(rates_df_component, fill_value=0)

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


class EquilibriumIonRateMatrix:
    def __init__(self):
        """Constructs ionization rate matrices from precomputed rate data."""

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
        cell : int
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
                    rate.index.get_level_values("ion_number_destination"),
                    rate.index.get_level_values("ion_number_source"),
                ),
            ),
            shape=(ion_states, ion_states),
        )

    def solve(
        self,
        photoion_rates_df: pd.DataFrame,
        recomb_rates_df: pd.DataFrame,
        collisional_ionization_rates_df: pd.DataFrame,
        collision_recombination_rates_df: pd.DataFrame,
        charge_conservation=False,
    ):
        """Compute the ionization rate matrix.

        Parameters
        ----------
        photoion_rates_df : pd.DataFrame
            Precomputed photoionization rates.
        recomb_rates_df : pd.DataFrame
            Precomputed radiative recombination rates.
        collisional_ionization_rates_df : pd.DataFrame
            Precomputed collisional ionization rates.
        collision_recombination_rates_df : pd.DataFrame
            Precomputed collisional recombination rates.
        charge_conservation : bool, optional
            Whether to include a charge conservation row in the rate matrix.

        Returns
        -------
        pd.DataFrame
            A DataFrame of rate matrices indexed by atomic number and ion number,
            with each column being a cell. Each entry is a numpy array.
        """
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
            index=list(grouped_photoion_rates_df.groups.keys()),
            columns=photoion_rates_df.columns,
        )

        for atomic_number in grouped_photoion_rates_df.groups.keys():
            photoion_rates = grouped_photoion_rates_df.get_group(atomic_number)
            recomb_rates = grouped_recomb_rates_df.get_group(atomic_number)
            coll_ion_rates = grouped_collisional_ionization_rates_df.get_group(
                atomic_number
            )
            recomb_ion_rates = (
                grouped_collisional_recombination_rates_df.get_group(
                    atomic_number
                )
            )
            max_ion_number = max(
                photoion_rates.index.get_level_values(
                    "ion_number_destination"
                ).max(),
                photoion_rates.index.get_level_values(
                    "ion_number_source"
                ).max(),
            )
            ion_states = int(max_ion_number) + 1
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
                matrix_array[1, :] = 1
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


class LTEIonRateMatrix:
    """Constructs ionization rate matrices based on LTE ratios."""

    def __init__(self):
        """Initialize the solver."""

    @staticmethod
    def _prepare_phi(phi, ion_index):
        """Prepare phi (Saha ratios) by reindexing to full ionization structure.

        Parameters
        ----------
        phi : pd.DataFrame
            Saha ratios indexed by (atomic_number, ion_number), columns are cells.
        ion_index : pd.MultiIndex
            Full ionization state index.

        Returns
        -------
        pd.DataFrame
            Phi reindexed to ion_index, with NaNs filled.
        """
        # Check for NaNs
        no_nans = pd.isnull(phi).sum().sum()
        if no_nans:
            phi = phi.fillna(phi.min().min())

        # Zero phi values cause numerical issues. Replace with small values.
        phi_min = phi[phi > 0.0].min().min()
        phi = phi.copy()
        phi[phi == 0.0] = 1.0e-10 * phi_min

        # Shift ion number by 1 and reindex to match ion population structure
        atomic_number = phi.index.get_level_values(0).values
        ion_number = phi.index.get_level_values(1).values
        new_index = pd.MultiIndex.from_arrays([atomic_number, ion_number - 1])
        phi_prep = phi.set_index(new_index).reindex(ion_index).fillna(0.0)
        return phi_prep

    @staticmethod
    def _get_number_conservation_index(ion_index):
        """Get indices for number conservation constraint row.

        Parameters
        ----------
        ion_index : pd.MultiIndex
            Index with (atomic_number, ion_number).

        Returns
        -------
        tuple
            (row_indices, col_indices) for setting number conservation row.
        """
        atomic_number = np.unique(
            ion_index.get_level_values(0).values.astype(int)
        )
        sum1 = (atomic_number + 1).cumsum() - 1
        index1 = np.concatenate(
            [np.ones(j + 1, dtype=int) * i for i, j in zip(sum1, atomic_number)]
        )
        index2 = np.arange(len(ion_index), dtype=int)
        return (index1, index2)

    def solve(
        self,
        phi,
        partition_function,
        electron_number_density,
        charge_conservation=False,
    ):
        """Compute ionization rate matrices from LTE Saha ratios.

        This produces a single rate matrix for each cell that governs ionization equilibrium
        for all species. The rate matrix includes electron density as a semi-variable,
        with number and charge conservation constraints.

        Parameters
        ----------
        phi : pd.DataFrame
            Saha ratios indexed by (atomic_number, ion_number), columns are cells.
        partition_function : pd.DataFrame
            Partition functions indexed by (atomic_number, ion_number), columns are cells.
        electron_number_density : np.ndarray
            Electron number density indexed by cell.
        charge_conservation : bool, optional
            Whether to include a charge conservation row in the rate matrix.

        Returns
        -------
        list of np.ndarray
            List of rate matrices, one per cell. Each matrix has shape
            (num_ions + 1, num_ions + 1) where num_ions is the total number of ion states
            and +1 accounts for electron density.
        """
        ion_index = partition_function.index

        # Prepare phi for use in rate matrix construction
        phi_prep = self._prepare_phi(phi, ion_index)

        # Get constraint information
        number_conservation_index = self._get_number_conservation_index(
            ion_index
        )

        # Construct rate matrix for each cell
        rate_matrices = pd.DataFrame(
            index=phi.index.get_level_values(0),
            columns=phi_prep.columns,
        )

        for atomic_number in phi.index.get_level_values(0):
            ion_states = atomic_number + 1
            for shell in range(len(phi_prep.columns)):
                # Get Saha ratios for this cell
                lte_diag = (
                    -phi_prep[shell].values / electron_number_density[shell]
                )
                lte_offdiag = (lte_diag != 0).astype(float)[:-1]

                matrix_array = np.diag(lte_diag) + np.diag(lte_offdiag, k=1)

                matrix_array[number_conservation_index] = 1.0
                if charge_conservation:
                    charge_conservation_row = np.hstack(
                        (np.arange(0, ion_states), -1)
                    )
                    matrix_array = np.pad(matrix_array, ((0, 0), (0, 1)))
                    matrix_array = np.vstack(
                        (charge_conservation_row, matrix_array)
                    )

                rate_matrices.loc[atomic_number, shell] = matrix_array

        return rate_matrices
