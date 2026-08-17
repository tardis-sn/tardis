import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.sparse import coo_matrix

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    EstimatedPhotoionizationRateSolver,
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_ion_population_to_level_population,
)
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
    PlanckianRadiationField,
)


class RateMatrix:
    """Build bound-bound rate matrices from rate solvers."""

    def __init__(
        self,
        radiative_rate_solver: RadiativeRatesSolver,
        electron_rate_solver: ThermalCollisionalRateSolver,
        levels: pd.DataFrame,
    ) -> None:
        """Construct a rate matrix from explicit bound-bound rate owners.

        Parameters
        ----------
        radiative_rate_solver : RadiativeRatesSolver
            Solver for radiative transition rates.
        electron_rate_solver : ThermalCollisionalRateSolver
            Solver for electron-dependent transition rates.
        levels : pd.DataFrame
            DataFrame of energy levels.
        """
        self.radiative_rate_solver = radiative_rate_solver
        self.electron_rate_solver = electron_rate_solver
        self.levels = levels

    def assemble_matrices(
        self,
        j_blues: pd.DataFrame,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        beta_sobolev: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """Assemble column-conserving bound-bound rate matrices.

        The returned matrices contain the rate equations, including
        their column-conserving diagonals, but do not contain the normalization
        row used by :meth:`solve`.  ``j_blues`` and ``beta_sobolev`` are used
        together for radiative transitions, so a residual evaluation can
        rebuild the matrix at a candidate Sobolev state.
        """
        rates_df_list = [
            self.radiative_rate_solver.solve(
                j_blues, beta_sobolev
            )
        ]
        rates_df_list.append(
            self.electron_rate_solver.solve(
                thermal_electron_energy_distribution.temperature
            )
        )

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
            rates_df * thermal_electron_energy_distribution.number_density.value
            if idx > 0
            else rates_df
            for idx, rates_df in enumerate(rates_df_list)
        ]

        rates_df = sum(rates_df_list)

        grouped_rates_df = rates_df.groupby(
            level=("atomic_number", "ion_number")
        )

        rate_matrices: dict[tuple[int, int], npt.NDArray[np.float64]] = {}

        for species_id, rates in grouped_rates_df:
            number_of_levels = self.levels.energy.loc[species_id].count()
            matrices = np.empty(
                (len(rates.columns), number_of_levels, number_of_levels)
            )
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
                matrices[shell] = matrix_array
            rate_matrices[species_id] = matrices

        rate_matrix_array = np.empty(
            (len(rate_matrices), len(rates_df.columns)), dtype=object
        )
        for species_idx, matrices in enumerate(rate_matrices.values()):
            rate_matrix_array[species_idx] = list(matrices)

        return pd.DataFrame(
            rate_matrix_array,
            index=pd.MultiIndex.from_tuples(
                rate_matrices, names=["atomic_number", "ion_number"]
            ),
            columns=pd.Index(rates_df.columns, dtype=object),
        )

    def solve(
        self,
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
    ) -> pd.DataFrame:
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
        j_blues = pd.DataFrame(
            radiation_field.calculate_mean_intensity(
                self.radiative_rate_solver.einstein_coefficients.nu.to_numpy()
            ),
            index=self.radiative_rate_solver.einstein_coefficients.index,
        )
        rate_matrices = self.assemble_matrices(
            j_blues,
            thermal_electron_energy_distribution,
        )
        normalized_matrices = rate_matrices.copy(deep=True)
        for species_id in normalized_matrices.index:
            for shell in normalized_matrices.columns:
                normalized_matrices.at[species_id, shell][0, :] = 1.0
        return normalized_matrices


class IonRateMatrix:
    """Build ionization rate matrices from radiative and collisional rates."""

    def __init__(
        self,
        radiative_ionization_rate_solver: AnalyticPhotoionizationRateSolver
        | EstimatedPhotoionizationRateSolver,
        collisional_ionization_rate_solver: CollisionalIonizationRateSolver,
    ) -> None:
        """Construct an ionization rate matrix solver.

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
        self.ion_population_index = pd.MultiIndex.from_tuples(
            [], names=["atomic_number", "ion_number"]
        )

    def __calculate_total_grouped_rates(
        self, rates_df: pd.DataFrame
    ) -> pd.core.groupby.DataFrameGroupBy:
        """Calculate grouped total ionization and recombination rates.

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

    def __construct_rate_matrix(
        self, rate: pd.DataFrame, cell: int, ion_states: int
    ) -> coo_matrix:
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
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        lte_level_population: pd.DataFrame,
        level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        level_to_continuum_saha_factor: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
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
        level_to_continuum_saha_factor : pandas.DataFrame, optional
            Density-independent Lucy level-to-continuum Saha factor. When
            omitted, retain the existing LTE-population-derived behavior.

        Returns
        -------
        pandas.DataFrame
            Rate matrices indexed by atomic number, with each column being a
            shell.
        """
        if isinstance(
            self.radiative_ionization_rate_solver,
            EstimatedPhotoionizationRateSolver,
        ):
            photoion_rates_df, recomb_rates_df = (
                self.radiative_ionization_rate_solver.solve(
                    thermal_electron_energy_distribution,
                    level_population,
                    ion_population,
                )
            )
        else:
            photoion_rates_df, recomb_rates_df = (
                self.radiative_ionization_rate_solver.solve(
                    radiation_field,
                    thermal_electron_energy_distribution,
                    lte_level_population,
                    level_population,
                    lte_ion_population,
                    ion_population,
                    partition_function,
                    boltzmann_factor,
                    level_to_continuum_saha_factor,
                )
            )

        if level_to_continuum_saha_factor is None:
            lte_ion_population = reindex_ion_population_to_level_population(
                lte_ion_population, lte_level_population
            )
            level_to_continuum_saha_factor = lte_level_population / (
                lte_ion_population.values
                * thermal_electron_energy_distribution.number_density.value
            )

        collisional_ionization_rates_df, collision_recombination_rates_df = (
            self.collisional_ionization_rate_solver.solve(
                thermal_electron_energy_distribution,
                level_to_continuum_saha_factor,
                partition_function,
                boltzmann_factor,
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

        rate_matrices = {}

        atomic_numbers = list(grouped_photoion_rates_df.groups.keys())
        self.ion_population_index = pd.MultiIndex.from_tuples(
            [
                (atomic_number, ion_number)
                for atomic_number in atomic_numbers
                for ion_number in range(atomic_number + 1)
            ],
            names=["atomic_number", "ion_number"],
        )

        for atomic_number in atomic_numbers:
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
            ion_states = atomic_number + 1
            matrix_arrays = np.empty(
                (len(photoion_rates.columns), ion_states, ion_states)
            )
            for shell_idx, shell in enumerate(photoion_rates.columns):
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
                matrix_arrays[shell_idx] = matrix_array
            rate_matrices[atomic_number] = matrix_arrays

        rate_matrix_array = np.empty(
            (len(rate_matrices), len(photoion_rates_df.columns)), dtype=object
        )
        for atomic_number_idx, matrices in enumerate(rate_matrices.values()):
            rate_matrix_array[atomic_number_idx] = list(matrices)

        return pd.DataFrame(
            rate_matrix_array,
            index=pd.Index(rate_matrices, name="atomic_number"),
            columns=photoion_rates_df.columns,
        )
