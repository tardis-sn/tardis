import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.inputs import BoundBoundMatrixRates
from tardis.plasma.equilibrium.matrix_assembly import (
    construct_rate_matrices,
    normalize_rate_matrices,
    sum_duplicate_rates,
    sum_rate_frames,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    EstimatedPhotoionizationRateSolver,
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.util import (
    align_ion_population_to_level_population,
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

    def prepare_species_rates(
        self,
        j_blues: pd.DataFrame,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        species: tuple[int, int],
        line_index: pd.Index,
    ) -> tuple[BoundBoundMatrixRates, ...]:
        """Prepare fixed transition arrays for every shell of one species."""
        radiative_rates = self.radiative_rate_solver.solve(j_blues)
        collisional_rates = self.electron_rate_solver.solve(
            thermal_electron_energy_distribution.temperature
        )
        radiative_species_mask = (
            radiative_rates.index.get_level_values("atomic_number")
            == species[0]
        ) & (
            radiative_rates.index.get_level_values("ion_number") == species[1]
        )
        collisional_species_mask = (
            collisional_rates.index.get_level_values("atomic_number")
            == species[0]
        ) & (
            collisional_rates.index.get_level_values("ion_number")
            == species[1]
        )
        radiative_rates = radiative_rates.loc[radiative_species_mask]
        collisional_rates = collisional_rates.loc[collisional_species_mask]
        transition_index = pd.MultiIndex.from_tuples(
            sorted(set(radiative_rates.index) | set(collisional_rates.index)),
            names=radiative_rates.index.names,
        )
        radiative_rate_coefficients = radiative_rates.reindex(
            transition_index, fill_value=0.0
        ).to_numpy(dtype=np.float64)
        collisional_rates = (
            collisional_rates.reindex(transition_index, fill_value=0.0)
            .set_axis(j_blues.columns, axis="columns")
            .to_numpy(dtype=np.float64)
            * thermal_electron_energy_distribution.number_density.value[
                None, :
            ]
        )

        source_level_idx = transition_index.get_level_values(
            "level_number_source"
        ).to_numpy(dtype=np.int64)
        destination_level_idx = transition_index.get_level_values(
            "level_number_destination"
        ).to_numpy(dtype=np.int64)
        lower_level_idx = np.minimum(source_level_idx, destination_level_idx)
        upper_level_idx = np.maximum(source_level_idx, destination_level_idx)
        radiative_line_index = pd.MultiIndex.from_arrays(
            [
                transition_index.get_level_values("atomic_number"),
                transition_index.get_level_values("ion_number"),
                lower_level_idx,
                upper_level_idx,
            ],
            names=line_index.names,
        )
        beta_line_idx = line_index.get_indexer(radiative_line_index).astype(
            np.int64
        )
        beta_line_idx[
            np.all(radiative_rate_coefficients == 0.0, axis=1)
        ] = -1
        number_of_levels = self.levels.energy.loc[species].count()
        return tuple(
            BoundBoundMatrixRates(
                number_of_levels,
                source_level_idx,
                destination_level_idx,
                radiative_rate_coefficients[:, shell_idx],
                collisional_rates[:, shell_idx],
                beta_line_idx,
            )
            for shell_idx in range(len(j_blues.columns))
        )

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
            self.radiative_rate_solver.solve(j_blues, beta_sobolev)
        ]
        rates_df_list.append(
            self.electron_rate_solver.solve(
                thermal_electron_energy_distribution.temperature
            )
        )

        rates_df = sum_rate_frames(
            rates_df_list,
            multipliers=[
                1.0,
                thermal_electron_energy_distribution.number_density.value,
            ],
        )

        grouped_rates_df = rates_df.groupby(
            level=("atomic_number", "ion_number")
        )

        rate_matrices: dict[tuple[int, int], npt.NDArray[np.float64]] = {}

        for species_id, rates in grouped_rates_df:
            number_of_levels = self.levels.energy.loc[species_id].count()
            matrices = construct_rate_matrices(
                rates,
                (number_of_levels, number_of_levels),
                "level_number_source",
                "level_number_destination",
            )
            for matrix in matrices:
                np.fill_diagonal(matrix, -np.sum(matrix, axis=0))
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


def assemble_ion_rate_matrices(
    photoion_rates_df: pd.DataFrame,
    recombination_rates_df: pd.DataFrame,
    collisional_ionization_rates_df: pd.DataFrame,
    collisional_recombination_rates_df: pd.DataFrame,
    lte_ionization_factor: pd.DataFrame | None = None,
    electron_density: npt.NDArray[np.float64] | None = None,
) -> tuple[pd.DataFrame, pd.MultiIndex]:
    """Assemble normalized ionization matrices from collected rates."""
    aggregated_rates = [
        sum_duplicate_rates(rates)
        for rates in (
            photoion_rates_df,
            recombination_rates_df,
            collisional_ionization_rates_df,
            collisional_recombination_rates_df,
        )
    ]
    rate_matrices: dict[int, npt.NDArray[np.float64]] = {}
    atomic_numbers = aggregated_rates[0].index.get_level_values(
        "atomic_number"
    ).unique()
    for atomic_number in atomic_numbers:
        rates = sum_rate_frames(
            [
                aggregated.xs(
                    atomic_number, level="atomic_number", drop_level=False
                )
                for aggregated in aggregated_rates
            ]
        )
        rate_matrices[atomic_number] = normalize_rate_matrices(
            construct_rate_matrices(
                rates,
                (atomic_number + 1, atomic_number + 1),
                "ion_number_source",
                "ion_number_destination",
            )
        )

    if lte_ionization_factor is not None:
        if electron_density is None:
            raise ValueError("electron_density is required with LTE rates")
        for atomic_number in lte_ionization_factor.index.get_level_values(
            "atomic_number"
        ).unique():
            if atomic_number in rate_matrices:
                continue
            ion_states = atomic_number + 1
            ionization_factor = lte_ionization_factor.loc[
                atomic_number, photoion_rates_df.columns
            ]
            matrix_arrays = np.zeros(
                (len(photoion_rates_df.columns), ion_states, ion_states)
            )
            for shell_idx, shell in enumerate(photoion_rates_df.columns):
                matrix_array = matrix_arrays[shell_idx]
                factors = ionization_factor[shell].to_numpy()
                ion_count = min(len(factors), ion_states - 1)
                source_ions = np.arange(ion_count)
                matrix_array[source_ions + 1, source_ions] = factors[:ion_count]
                matrix_array[source_ions, source_ions + 1] = electron_density[
                    shell_idx
                ]
            rate_matrices[atomic_number] = normalize_rate_matrices(matrix_arrays)

    atomic_numbers = sorted(rate_matrices)
    ion_population_index = pd.MultiIndex.from_tuples(
        [
            (atomic_number, ion_number)
            for atomic_number in atomic_numbers
            for ion_number in range(atomic_number + 1)
        ],
        names=["atomic_number", "ion_number"],
    )
    rate_matrix_array = np.empty(
        (len(rate_matrices), len(photoion_rates_df.columns)), dtype=object
    )
    for atomic_number_idx, atomic_number in enumerate(atomic_numbers):
        rate_matrix_array[atomic_number_idx] = list(rate_matrices[atomic_number])
    return (
        pd.DataFrame(
            rate_matrix_array,
            index=pd.Index(atomic_numbers, name="atomic_number"),
            columns=photoion_rates_df.columns,
        ),
        ion_population_index,
    )


class AnalyticIonRateMatrix:
    """Build ionization matrices from analytic radiative rates."""

    def __init__(
        self,
        radiative_ionization_rate_solver: AnalyticPhotoionizationRateSolver,
        collisional_ionization_rate_solver: CollisionalIonizationRateSolver,
    ) -> None:
        self.radiative_ionization_rate_solver = radiative_ionization_rate_solver
        self.collisional_ionization_rate_solver = collisional_ionization_rate_solver
        self.ion_population_index = pd.MultiIndex.from_tuples(
            [], names=["atomic_number", "ion_number"]
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
        lte_ionization_factor: pd.DataFrame | None = None,
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
        if level_to_continuum_saha_factor is None:
            lte_ion_population = align_ion_population_to_level_population(
                lte_ion_population, lte_level_population
            )
            level_to_continuum_saha_factor = lte_level_population / (
                lte_ion_population.values
                * thermal_electron_energy_distribution.number_density.value
            )

        photoion_rates_df, recomb_rates_df = self.radiative_ionization_rate_solver.solve(
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
        collisional_ionization_rates_df, collision_recombination_rates_df = (
            self.collisional_ionization_rate_solver.solve(
                thermal_electron_energy_distribution,
                level_to_continuum_saha_factor,
                partition_function,
                boltzmann_factor,
            )
        )
        rate_matrices, self.ion_population_index = assemble_ion_rate_matrices(
            photoion_rates_df,
            recomb_rates_df,
            collisional_ionization_rates_df,
            collision_recombination_rates_df,
            lte_ionization_factor=lte_ionization_factor,
            electron_density=thermal_electron_energy_distribution.number_density.to_value(
                "cm^-3"
            ),
        )
        return rate_matrices


class EstimatedIonRateMatrix:
    """Build ionization matrices from fixed Monte Carlo estimator rates."""

    def __init__(
        self,
        radiative_ionization_rate_solver: EstimatedPhotoionizationRateSolver,
        collisional_ionization_rate_solver: CollisionalIonizationRateSolver,
        lte_ionization_factor: pd.DataFrame | None = None,
    ) -> None:
        self.radiative_ionization_rate_solver = radiative_ionization_rate_solver
        self.collisional_ionization_rate_solver = collisional_ionization_rate_solver
        self.lte_ionization_factor = lte_ionization_factor
        self.ion_population_index = pd.MultiIndex.from_tuples(
            [], names=["atomic_number", "ion_number"]
        )
        self.prepared_matrix_coefficients: dict[
            int, tuple[npt.NDArray[np.float64], ...]
        ] | None = None
        self.prepared_columns: pd.Index | None = None

    def prepare(
        self,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        level_population: pd.DataFrame,
        ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        level_to_continuum_saha_factor: pd.DataFrame,
        lte_ionization_factor: pd.DataFrame | None = None,
    ) -> None:
        """Prepare density-polynomial ion matrices for charge roots."""
        unit_density_distribution = ThermalElectronEnergyDistribution(
            thermal_electron_energy_distribution.energy,
            thermal_electron_energy_distribution.temperature,
            np.ones_like(
                thermal_electron_energy_distribution.number_density.value
            )
            * u.cm**-3,
        )
        photoionization_rate, recombination_rate = (
            self.radiative_ionization_rate_solver.solve(
                unit_density_distribution,
                level_population,
                ion_population,
                level_to_continuum_saha_factor,
            )
        )
        collisional_ionization_rate, collisional_recombination_rate = (
            self.collisional_ionization_rate_solver.solve(
                unit_density_distribution,
                level_to_continuum_saha_factor,
                partition_function,
                boltzmann_factor,
                level_population=level_population,
                ion_population=ion_population,
            )
        )
        aggregated_rates = [
            sum_duplicate_rates(rates)
            for rates in (
                photoionization_rate,
                recombination_rate,
                collisional_ionization_rate,
                collisional_recombination_rate,
            )
        ]
        matrix_coefficients: dict[
            int, tuple[npt.NDArray[np.float64], ...]
        ] = {}
        atomic_numbers = aggregated_rates[0].index.get_level_values(
            "atomic_number"
        ).unique()
        for atomic_number in atomic_numbers:
            constant_rates = aggregated_rates[0].xs(
                atomic_number, level="atomic_number", drop_level=False
            )
            linear_rates = sum_rate_frames(
                [
                    rates.xs(
                        atomic_number,
                        level="atomic_number",
                        drop_level=False,
                    )
                    for rates in aggregated_rates[1:3]
                ]
            )
            quadratic_rates = aggregated_rates[3].xs(
                atomic_number, level="atomic_number", drop_level=False
            )
            coefficients = []
            for rates in (constant_rates, linear_rates, quadratic_rates):
                matrices = construct_rate_matrices(
                    rates,
                    (atomic_number + 1, atomic_number + 1),
                    "ion_number_source",
                    "ion_number_destination",
                )
                for matrix in matrices:
                    np.fill_diagonal(matrix, -np.sum(matrix, axis=0))
                coefficients.append(matrices)
            coefficients[0][:, 1, :] = 1.0
            coefficients[1][:, 1, :] = 0.0
            coefficients[2][:, 1, :] = 0.0
            matrix_coefficients[atomic_number] = tuple(coefficients)

        effective_lte_ionization_factor = (
            lte_ionization_factor
            if lte_ionization_factor is not None
            else self.lte_ionization_factor
        )
        if effective_lte_ionization_factor is not None:
            for atomic_number in (
                effective_lte_ionization_factor.index.get_level_values(
                    "atomic_number"
                ).unique()
            ):
                if atomic_number in matrix_coefficients:
                    continue
                ion_states = atomic_number + 1
                constant_matrices = np.zeros(
                    (len(level_population.columns), ion_states, ion_states)
                )
                linear_matrices = np.zeros_like(constant_matrices)
                ionization_factor = effective_lte_ionization_factor.loc[
                    atomic_number, level_population.columns
                ]
                for shell_idx, shell in enumerate(level_population.columns):
                    factors = ionization_factor[shell].to_numpy()
                    ion_count = min(len(factors), ion_states - 1)
                    source_ions = np.arange(ion_count)
                    constant_matrices[
                        shell_idx, source_ions + 1, source_ions
                    ] = factors[:ion_count]
                    linear_matrices[
                        shell_idx, source_ions, source_ions + 1
                    ] = 1.0
                for matrices in (constant_matrices, linear_matrices):
                    for matrix in matrices:
                        np.fill_diagonal(matrix, -np.sum(matrix, axis=0))
                constant_matrices[:, 1, :] = 1.0
                linear_matrices[:, 1, :] = 0.0
                matrix_coefficients[atomic_number] = (
                    constant_matrices,
                    linear_matrices,
                    np.zeros_like(constant_matrices),
                )

        atomic_numbers = sorted(matrix_coefficients)
        self.ion_population_index = pd.MultiIndex.from_tuples(
            [
                (atomic_number, ion_number)
                for atomic_number in atomic_numbers
                for ion_number in range(atomic_number + 1)
            ],
            names=["atomic_number", "ion_number"],
        )
        self.prepared_matrix_coefficients = matrix_coefficients
        self.prepared_columns = level_population.columns.copy()

    def solve_prepared(
        self,
        electron_density: npt.NDArray[np.float64],
        columns: pd.Index,
    ) -> pd.DataFrame:
        """Evaluate prepared ion matrices at trial electron densities."""
        if (
            self.prepared_matrix_coefficients is None
            or self.prepared_columns is None
        ):
            raise RuntimeError("Ion rate coefficients have not been prepared.")
        shell_indices = self.prepared_columns.get_indexer(columns)
        if np.any(shell_indices < 0):
            raise ValueError("Prepared ion matrices do not contain every shell.")
        density = np.asarray(electron_density, dtype=np.float64)[:, None, None]
        atomic_numbers = sorted(self.prepared_matrix_coefficients)
        rate_matrix_array = np.empty(
            (len(atomic_numbers), len(columns)), dtype=object
        )
        for atomic_number_idx, atomic_number in enumerate(atomic_numbers):
            constant, linear, quadratic = self.prepared_matrix_coefficients[
                atomic_number
            ]
            matrices = (
                constant[shell_indices]
                + density * linear[shell_indices]
                + density**2 * quadratic[shell_indices]
            )
            rate_matrix_array[atomic_number_idx] = list(matrices)
        return pd.DataFrame(
            rate_matrix_array,
            index=pd.Index(atomic_numbers, name="atomic_number"),
            columns=columns,
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
        level_to_continuum_saha_factor: pd.DataFrame,
        lte_ionization_factor: pd.DataFrame | None = None,
    ) -> pd.DataFrame:
        """Compute the ionization rate matrix from fixed estimators."""
        del radiation_field, lte_level_population, lte_ion_population
        self.prepare(
            thermal_electron_energy_distribution,
            level_population,
            ion_population,
            partition_function,
            boltzmann_factor,
            level_to_continuum_saha_factor,
            lte_ionization_factor,
        )
        return self.solve_prepared(
            thermal_electron_energy_distribution.number_density.to_value(
                "cm^-3"
            ),
            level_population.columns,
        )
