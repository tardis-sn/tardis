import logging

import astropy.units as u
import numpy as np
import numpy.typing as npt
import pandas as pd

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.charge_conservation import (
    ChargeConservationSolver,
)
from tardis.plasma.equilibrium.population_state import (
    PopulationState,
    SingleElementPopulationState,
)
from tardis.plasma.equilibrium.rate_matrix import (
    IonLevelRateMatrixSet,
)
from tardis.plasma.exceptions import PlasmaIonizationError
from tardis.plasma.properties.level_population import LevelNumberDensity

logger = logging.getLogger(__name__)

LOWER_ION_LEVEL_H = 0
_POPULATION_TOLERANCE = 1e-10


def get_lower_ion_level_index(
    level_population: pd.DataFrame,
) -> npt.NDArray[np.bool_]:
    """Return the mask for levels belonging to the lower ionization stage."""
    return (
        level_population.index.get_level_values("ion_number")
        == LOWER_ION_LEVEL_H
    )


def get_upper_ion_population_index(
    ion_population: pd.DataFrame,
) -> npt.NDArray[np.bool_]:
    """Return the mask for ion populations above the lower ionization stage."""
    return (
        ion_population.index.get_level_values("ion_number") > LOWER_ION_LEVEL_H
    )


def calculate_level_to_ion_population_factor(
    level_population_at_lower_ion: pd.DataFrame,
    ion_population_at_upper_ion: pd.DataFrame,
    electron_number_density: u.Quantity,
) -> pd.DataFrame:
    """Compute the Lucy 2003 Eq. 14 level-to-ion population factor."""
    density = electron_number_density.value
    # The factor is multiplied by one or two powers of density downstream.
    # Use its finite density-independent ratio at the zero-density endpoint.
    density_for_factor = np.where(density == 0.0, 1.0, density)
    return level_population_at_lower_ion / (
        ion_population_at_upper_ion.values * density_for_factor
    )


class AnalyticIonPopulationSolver:
    """Solve ion populations using analytic radiative rates."""

    def __init__(
        self,
        rate_matrix_solver: object,
        photoionization_rate_solver: object,
        collisional_ionization_rate_solver: object,
        elemental_number_density: pd.DataFrame,
        *,
        radiation_field: object | None = None,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution
        | None = None,
        lte_level_population: pd.DataFrame | None = None,
        estimated_level_population: pd.DataFrame | None = None,
        lte_ion_population: pd.DataFrame | None = None,
        estimated_ion_population: pd.DataFrame | None = None,
        partition_function: pd.DataFrame | float | None = None,
        boltzmann_factor: pd.DataFrame | None = None,
        max_solver_iterations: int = 100,
        tolerance: float = 1e-14,
    ) -> None:
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver : EquilibriumIonRateMatrix
        photoionization_rate_solver : RadiativeIonizationRateSolver
        collisional_ionization_rate_solver : CollisionalIonizationRateSolver
        elemental_number_density : pandas.DataFrame
            Elemental number density indexed by atomic number and columned by
            shell.
        radiation_field : object, optional
            Radiation field fixed across population evaluations.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution, optional
            Electron state used by the population solve.
        lte_level_population, estimated_level_population : pandas.DataFrame, optional
            LTE and estimated level populations.
        lte_ion_population, estimated_ion_population : pandas.DataFrame, optional
            LTE and estimated ion populations.
        partition_function, boltzmann_factor : pandas.DataFrame, optional
            Partition functions and level Boltzmann factors.
        max_solver_iterations : int, optional
            Maximum number of iterations for the ion population solver, by default 100.
        tolerance : float, optional
            Convergence tolerance for ion and electron populations, by default
            ``1e-14``.
        """
        self.rate_matrix_solver = rate_matrix_solver
        self.photoionization_rate_solver = photoionization_rate_solver
        self.collisional_ionization_rate_solver = (
            collisional_ionization_rate_solver
        )
        self.elemental_number_density = elemental_number_density
        self.radiation_field = radiation_field
        self.thermal_electron_energy_distribution = (
            thermal_electron_energy_distribution
        )
        self.lte_level_population = lte_level_population
        self.estimated_level_population = estimated_level_population
        self.lte_ion_population = lte_ion_population
        self.estimated_ion_population = estimated_ion_population
        self.partition_function = partition_function
        self.boltzmann_factor = boltzmann_factor
        self.max_solver_iterations = max_solver_iterations
        self.tolerance = tolerance

    @staticmethod
    def calculate_balance_vector(
        elemental_number_density_in_cell,
        rate_matrix_index,
        charge_conservation=False,
    ):
        """Construct the balance vector for the NLTE ionization equations.

        Combines all solution-vector blocks.

        Parameters
        ----------
        elemental_number_density : pandas.Series
            Elemental number densities indexed by atomic number.
        rate_matrix_index : pandas.MultiIndex
            (atomic_number, ion_number)
        charge_conservation : bool, optional
            Whether to include a charge conservation row.

        Returns
        -------
        numpy.array
            Solution vector for the NLTE ionization solver.
        """
        balance_array = [
            np.append(
                np.zeros(
                    (
                        rate_matrix_index.get_level_values("atomic_number")
                        == atomic_number
                    ).sum()
                ),
                elemental_number_density_in_cell.loc[atomic_number],
            )
            for atomic_number in elemental_number_density_in_cell.index
        ]

        if charge_conservation:
            balance_array.append(np.array([0.0]))

        return np.hstack(balance_array)

    @staticmethod
    def delta_quantity(estimated, solution):
        """Return the relative change between estimated and solved values."""
        return (estimated - solution) / solution

    @staticmethod
    def _build_elemental_population_solution(
        elemental_rate_matrices: IonLevelRateMatrixSet,
        elemental_number_density: pd.Series,
    ) -> SingleElementPopulationState:
        """Solve and unpack all shells in an elemental rate-matrix set.

        The normalized matrix is solved once per shell. The returned state
        vector is then unpacked exclusively through ``state_index``: explicit
        level positions become level populations, while all level positions
        belonging to one ion are summed to recover its total population. Ion
        stages without explicit levels are read from their single total-state
        position, which preserves the terminal bare-nucleus state.

        Absolute populations are reconstructed only after the normalized
        state is complete by multiplying each shell by the elemental number
        density. This keeps the linear solve independent of the abundance
        scale.

        Parameters
        ----------
        elemental_rate_matrices : IonLevelRateMatrixSet
            Element-normalized matrices and explicit state-index metadata from
            :meth:`~tardis.plasma.equilibrium.rate_matrix.IonRateMatrix.solve_ion_and_level`.
        elemental_number_density : pandas.Series
            Absolute number density of the element indexed by matrix shell.

        Returns
        -------
        SingleElementPopulationState
            Element-normalized level and ion fractions and reconstructed
            absolute populations.
        """
        normalized_rate_matrices = (
            elemental_rate_matrices.normalized_rate_matrices
        )
        state_index = elemental_rate_matrices.state_index
        atomic_number = state_index.atomic_number
        columns = normalized_rate_matrices.columns
        density = elemental_number_density.reindex(columns).astype(float)

        state_populations = pd.DataFrame(
            index=state_index.states,
            columns=columns,
            dtype=float,
        )
        for cell in columns:
            # The matrix state index is the single source of truth for the
            # shell vector; level and ion offsets must not be reconstructed.
            normalized_rate_matrix = normalized_rate_matrices.loc[
                atomic_number, cell
            ]
            normalization = np.zeros(normalized_rate_matrix.shape[0])
            normalization[0] = 1.0
            try:
                population = np.linalg.solve(
                    normalized_rate_matrix,
                    normalization,
                )
            except np.linalg.LinAlgError as exc:
                raise PlasmaIonizationError(
                    "Singular population matrix for atomic number "
                    f"{atomic_number}, shell {cell}"
                ) from exc

            raw_rate_matrix = (
                elemental_rate_matrices.raw_elemental_rate_matrices.loc[
                    atomic_number, cell
                ]
            )
            residual = raw_rate_matrix @ population
            row_scale = np.max(np.abs(raw_rate_matrix), axis=1)
            nonzero_rows = row_scale > 0.0
            scaled_residual = np.zeros_like(residual)
            scaled_residual[nonzero_rows] = (
                np.abs(residual[nonzero_rows]) / row_scale[nonzero_rows]
            )
            minimum_population = float(np.min(population))
            worst_row = int(np.argmax(scaled_residual))
            if (
                not np.isfinite(population).all()
                or minimum_population < -_POPULATION_TOLERANCE
                or not np.isclose(
                    population.sum(),
                    1.0,
                    rtol=0.0,
                    atol=_POPULATION_TOLERANCE,
                )
                or scaled_residual[worst_row] > _POPULATION_TOLERANCE
            ):
                raise PlasmaIonizationError(
                    "Invalid population solution for atomic number "
                    f"{atomic_number}, shell {cell}, row {worst_row}: "
                    f"minimum={minimum_population}, "
                    f"normalization={population.sum()}, "
                    f"scaled_residual={scaled_residual[worst_row]}"
                )
            population[population < 0.0] = 0.0
            population /= population.sum()
            state_populations[cell] = population

        level_positions = sorted(
            state_index.level_positions.items(), key=lambda item: item[1]
        )
        level_index = pd.MultiIndex.from_tuples(
            [
                (atomic_number, ion_number, level_number)
                for (ion_number, level_number), _ in level_positions
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        normalized_level_populations = pd.DataFrame(
            [
                state_populations.iloc[position].to_numpy()
                for _, position in level_positions
            ],
            index=level_index,
            columns=columns,
        )

        # Sum explicit level states to recover the total population of each
        # level-bearing ion. Bare and level-free ions already have one total
        # state and therefore require no summation.
        ion_numbers = sorted(
            set(state_index.ion_positions)
            | {ion_number for ion_number, _ in state_index.level_positions}
        )
        normalized_ion_populations = pd.DataFrame(
            index=pd.MultiIndex.from_tuples(
                [(atomic_number, ion_number) for ion_number in ion_numbers],
                names=["atomic_number", "ion_number"],
            ),
            columns=columns,
            dtype=float,
        )
        for ion_number in ion_numbers:
            ion_level_positions = [
                position
                for (
                    level_ion_number,
                    _,
                ), position in state_index.level_positions.items()
                if level_ion_number == ion_number
            ]
            if ion_level_positions:
                ion_population = state_populations.iloc[
                    ion_level_positions
                ].sum(axis=0)
            else:
                ion_population = state_populations.iloc[
                    state_index.ion_positions[ion_number]
                ]
            normalized_ion_populations.loc[(atomic_number, ion_number)] = (
                ion_population.to_numpy()
            )

        level_populations = normalized_level_populations.multiply(
            density, axis="columns"
        )
        ion_populations = normalized_ion_populations.multiply(
            density, axis="columns"
        )
        return SingleElementPopulationState(
            normalized_state_populations=state_populations,
            normalized_level_populations=normalized_level_populations,
            normalized_ion_populations=normalized_ion_populations,
            level_populations=level_populations,
            ion_populations=ion_populations,
            state_index=state_index,
        )

    def _calculate_rates(
        self,
        electron_energy_distribution: ThermalElectronEnergyDistribution,
        radiation_field: object,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
    ) -> tuple[pd.DataFrame, ...]:
        """Calculate analytic radiative and collisional rate frames."""
        lower_ion_level_index = lte_level_population.index[
            get_lower_ion_level_index(lte_level_population)
        ]
        upper_ion_population_index = lte_ion_population.index[
            get_upper_ion_population_index(lte_ion_population)
        ]
        photoion_rates_df, recomb_rates_df = (
            self.photoionization_rate_solver.solve(
                radiation_field,
                electron_energy_distribution,
                lte_level_population.loc[lower_ion_level_index],
                estimated_level_population.loc[lower_ion_level_index],
                lte_ion_population.loc[upper_ion_population_index],
                estimated_ion_population.loc[upper_ion_population_index],
                partition_function,
                boltzmann_factor,
            )
        )
        level_to_ion_population_factor = (
            calculate_level_to_ion_population_factor(
                lte_level_population.loc[lower_ion_level_index],
                lte_ion_population.loc[upper_ion_population_index],
                electron_energy_distribution.number_density,
            )
        )
        collisional_ionization_rates_df, collision_recombination_rates_df = (
            self.collisional_ionization_rate_solver.solve(
                electron_energy_distribution,
                level_to_ion_population_factor,
                partition_function,
                boltzmann_factor,
            )
        )
        return (
            photoion_rates_df,
            recomb_rates_df,
            collisional_ionization_rates_df,
            collision_recombination_rates_df,
        )

    def solve(self) -> tuple[pd.DataFrame, pd.Series]:
        """Solves the normalized ion population values from the rate matrices.

        Returns
        -------
        pd.DataFrame
            Normalized ion population values indexed by atomic number, ion
            number. Columns are cells.
        pd.Series
            Electron number density values. Index is cells.
        """
        new_electron_energy_distribution = (
            self.thermal_electron_energy_distribution
        )
        estimated_ion_population = self.estimated_ion_population
        for iteration in range(self.max_solver_iterations):
            (
                photoion_rates_df,
                recomb_rates_df,
                collisional_ionization_rates_df,
                collision_recombination_rates_df,
            ) = self._calculate_rates(
                new_electron_energy_distribution,
                self.radiation_field,
                self.lte_level_population,
                self.estimated_level_population,
                self.lte_ion_population,
                estimated_ion_population,
                self.partition_function,
                self.boltzmann_factor,
            )
            rates_matrices = self.rate_matrix_solver.solve(
                photoion_rates_df,
                recomb_rates_df,
                collisional_ionization_rates_df,
                collision_recombination_rates_df,
            )
            self.rates_matrices = rates_matrices
            solved_matrices = pd.DataFrame(
                index=rates_matrices.index,
                columns=rates_matrices.columns,
            )
            for cell in self.elemental_number_density.columns:
                balance_vector = self.calculate_balance_vector(
                    self.elemental_number_density[cell],
                    rates_matrices.index,
                )
                solved_matrix = rates_matrices[cell].apply(
                    np.linalg.solve, args=(balance_vector,)
                )
                solved_matrices[cell] = solved_matrix

            ion_population_solution = pd.DataFrame(
                np.vstack(solved_matrices.values[0]).T,
                index=estimated_ion_population.index,
                columns=rates_matrices.columns,
            )

            if (ion_population_solution < 0).any().any():
                ion_population_solution[ion_population_solution < 0] = 0.0

            electron_population_solution = (
                ion_population_solution
                * ion_population_solution.index.get_level_values("ion_number")
                .values[np.newaxis]
                .T
            ).sum()

            delta_ion = self.delta_quantity(
                estimated_ion_population, ion_population_solution
            )
            delta_electron = self.delta_quantity(
                new_electron_energy_distribution.number_density.value,
                electron_population_solution,
            )

            if (
                np.all(np.abs(delta_ion) < self.tolerance).any().any()
                and (np.abs(delta_electron) < self.tolerance).any().any()
            ):
                logger.info(
                    "Ion population solver converged after %d iterations.",
                    iteration + 1,
                )
                break

            estimated_ion_population = ion_population_solution
            new_electron_energy_distribution.number_density = (
                electron_population_solution.values * u.cm**-3
            )
        else:
            logger.warning(
                "Ion population solver did not converge after %d iterations.",
                iteration,
            )

        return ion_population_solution, electron_population_solution


class ChargeConservingAnalyticIonPopulationSolver(AnalyticIonPopulationSolver):
    """Solve analytic ion populations with global charge conservation."""

    def _solve_elemental_populations(
        self,
        rate_frames: tuple[pd.DataFrame, ...],
        columns: pd.Index,
    ) -> dict[int, SingleElementPopulationState]:
        """Solve each element covered by the supplied rate frames."""
        elemental_populations: dict[int, SingleElementPopulationState] = {}
        for atomic_number in self.elemental_number_density.index:
            elemental_rate_frames = tuple(
                rate_frame.loc[
                    rate_frame.index.get_level_values("atomic_number")
                    == atomic_number
                ]
                for rate_frame in rate_frames
            )
            if not any(
                not rate_frame.empty for rate_frame in elemental_rate_frames
            ):
                raise ValueError(
                    "Fixed-density rates do not cover atomic number "
                    f"{atomic_number}"
                )
            matrix_set = self.rate_matrix_solver.solve_ion_and_level(
                int(atomic_number),
                photoion_rates_df=elemental_rate_frames[0],
                recomb_rates_df=elemental_rate_frames[1],
                collisional_ionization_rates_df=elemental_rate_frames[2],
                collision_recombination_rates_df=elemental_rate_frames[3],
            )
            elemental_populations[int(atomic_number)] = (
                self._build_elemental_population_solution(
                    matrix_set,
                    self.elemental_number_density.loc[atomic_number, columns],
                )
            )
        return elemental_populations

    def solve_charge_state_at_electron_density(
        self,
        electron_number_density: pd.Series,
    ) -> PopulationState:
        """Solve the analytic state at a trial electron density."""
        columns = electron_number_density.index
        shell_positions = self.elemental_number_density.columns.get_indexer(
            columns
        )
        electron_distribution = ThermalElectronEnergyDistribution(
            self.thermal_electron_energy_distribution.energy,
            self.thermal_electron_energy_distribution.temperature[
                shell_positions
            ],
            electron_number_density.to_numpy() * u.cm**-3,
        )
        rate_frames = self._calculate_rates(
            electron_distribution,
            self.radiation_field,
            self.lte_level_population.loc[:, columns],
            self.estimated_level_population.loc[:, columns],
            self.lte_ion_population.loc[:, columns],
            self.estimated_ion_population.loc[:, columns],
            self.partition_function.loc[:, columns],
            self.boltzmann_factor.loc[:, columns],
        )
        rate_frames = tuple(
            rate_frame.reindex(columns=columns) for rate_frame in rate_frames
        )
        elemental_populations = self._solve_elemental_populations(
            rate_frames, columns
        )
        ion_number_density = pd.concat(
            [
                solution.ion_populations
                for solution in elemental_populations.values()
            ]
        ).sort_index()
        level_boltzmann_factor = self.boltzmann_factor.loc[:, columns].copy(
            deep=True
        )
        level_number_density = LevelNumberDensity(None).calculate(
            level_boltzmann_factor,
            ion_number_density,
            level_boltzmann_factor.index,
            self.partition_function.loc[:, columns],
        )
        level_number_density.columns = columns
        for solution in elemental_populations.values():
            if not solution.level_populations.empty:
                level_number_density.loc[solution.level_populations.index] = (
                    solution.level_populations
                )
        return PopulationState(
            electron_densities=electron_number_density.copy(),
            elemental_populations=elemental_populations,
            ion_number_density=ion_number_density,
            level_number_density=level_number_density,
            level_boltzmann_factor=level_boltzmann_factor,
        )

    def solve(self) -> tuple[pd.DataFrame, pd.Series]:
        """Return analytic ion populations at charge-conserving densities."""
        population_state = ChargeConservationSolver(
            self.elemental_number_density,
            self,
            max_solver_iterations=self.max_solver_iterations,
        ).solve()
        return (
            population_state.ion_number_density,
            population_state.electron_densities,
        )


class EstimatedIonPopulationSolver(AnalyticIonPopulationSolver):
    """Solve ion populations using Monte Carlo rate estimators."""

    def _calculate_rates(
        self,
        electron_energy_distribution: ThermalElectronEnergyDistribution,
        radiation_field: object,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
    ) -> tuple[pd.DataFrame, ...]:
        """Calculate estimator-backed radiative and collisional rate frames."""
        lower_ion_level_index = lte_level_population.index[
            get_lower_ion_level_index(lte_level_population)
        ]
        upper_ion_population_index = lte_ion_population.index[
            get_upper_ion_population_index(lte_ion_population)
        ]
        photoion_rates_df, recomb_rates_df = (
            self.photoionization_rate_solver.solve(
                electron_energy_distribution,
                radiation_field,
                estimated_level_population.loc[lower_ion_level_index],
                lte_level_population.loc[lower_ion_level_index],
                estimated_ion_population.loc[upper_ion_population_index],
                lte_ion_population.loc[upper_ion_population_index],
            )
        )
        level_to_ion_population_factor = (
            calculate_level_to_ion_population_factor(
                lte_level_population.loc[lower_ion_level_index],
                lte_ion_population.loc[upper_ion_population_index],
                electron_energy_distribution.number_density,
            )
        )
        collisional_ionization_rates_df, collision_recombination_rates_df = (
            self.collisional_ionization_rate_solver.solve(
                electron_energy_distribution,
                level_to_ion_population_factor,
                partition_function,
                boltzmann_factor,
            )
        )
        return (
            photoion_rates_df,
            recomb_rates_df,
            collisional_ionization_rates_df,
            collision_recombination_rates_df,
        )


class LTEIonPopulationSolver(AnalyticIonPopulationSolver):
    """Solve ion populations from LTE rate matrices."""

    def __init__(
        self,
        rate_matrix_solver: object,
        elemental_number_density: pd.DataFrame,
        max_solver_iterations: int = 100,
        charge_conservation: bool = False,
        tolerance: float = 1e-14,
    ) -> None:
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver : LTEIonRateMatrix
        elemental_number_density : pd.DataFrame
            Elemental number density. Index is atomic number, columns are cells.
        max_solver_iterations : int, optional
            Maximum number of iterations for the ion population solver.
        charge_conservation : bool, optional
            Whether to include charge conservation in the rate matrix.
        tolerance : float, optional
            Convergence tolerance for ion and electron populations.
        """
        super().__init__(
            rate_matrix_solver,
            None,
            None,
            elemental_number_density,
            max_solver_iterations=max_solver_iterations,
            tolerance=tolerance,
        )
        self.charge_conservation = charge_conservation

    def solve(
        self,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        lte_ion_population: pd.DataFrame,
        phi: pd.DataFrame,
        partition_function: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Solves the normalized ion population values from the rate matrices.

        Parameters
        ----------
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        lte_ion_population : pd.DataFrame
            LTE ion number density. Columns are cells.
        phi : pd.DataFrame
            Radiation field mean intensity. Columns are cells.
        partition_function : pd.DataFrame
            Partition function values. Columns are cells.

        Returns
        -------
        pd.DataFrame
            Normalized ion population values indexed by atomic number, ion
            number. Columns are cells.
        pd.Series
            Electron number density values. Index is cells.
        """
        new_electron_energy_distribution = thermal_electron_energy_distribution

        for iteration in range(self.max_solver_iterations):
            rates_matrices = self.rate_matrix_solver.solve(
                phi,
                partition_function,
                new_electron_energy_distribution.number_density.value,
                self.charge_conservation,
            )
            self.rates_matrices = rates_matrices
            solved_matrices = pd.DataFrame(
                index=rates_matrices.index,
                columns=rates_matrices.columns,
            )
            for cell in self.elemental_number_density.columns:
                balance_vector = self.calculate_balance_vector(
                    self.elemental_number_density[cell],
                    rates_matrices.index,
                    self.charge_conservation,
                )
                solved_matrix = rates_matrices[cell].apply(
                    np.linalg.solve, args=(balance_vector,)
                )
                solved_matrices[cell] = solved_matrix

            ion_population_solution = pd.DataFrame(
                np.vstack(solved_matrices.values[0]).T,
                index=lte_ion_population.index,
                columns=rates_matrices.columns,
            )

            if (ion_population_solution < 0).any().any():
                ion_population_solution[ion_population_solution < 0] = 0.0

            electron_population_solution = (
                ion_population_solution
                * ion_population_solution.index.get_level_values("ion_number")
                .values[np.newaxis]
                .T
            ).sum()

            delta_ion = self.delta_quantity(
                lte_ion_population, ion_population_solution
            )
            delta_electron = self.delta_quantity(
                new_electron_energy_distribution.number_density.value,
                electron_population_solution,
            )

            if (
                np.all(np.abs(delta_ion) < self.tolerance).any().any()
                and (np.abs(delta_electron) < self.tolerance).any().any()
            ):
                logger.info(
                    "Ion population solver converged after %d iterations.",
                    iteration + 1,
                )
                break

            lte_ion_population = ion_population_solution
            new_electron_energy_distribution.number_density = (
                electron_population_solution.values * u.cm**-3
            )
        else:
            logger.warning(
                "Ion population solver did not converge after %d iterations.",
                iteration,
            )

        return ion_population_solution, electron_population_solution
