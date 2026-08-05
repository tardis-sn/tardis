import logging
from collections.abc import Callable
from dataclasses import dataclass
from functools import partial

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
from tardis.plasma.equilibrium.rate_matrix import (
    ElementalStateIndex,
    IonLevelRateMatrixSet,
)

logger = logging.getLogger(__name__)

LOWER_ION_LEVEL_H = 0


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


@dataclass(frozen=True)
class ElementalPopulationSolution:
    """Element-normalized and absolute populations for one element.

    Attributes
    ----------
    normalized_state_populations : pandas.DataFrame
        Population for every explicit level and total-ion state in the
        matrix ordering described by ``state_index.states``. Columns are
        shells.
    normalized_level_populations : pandas.DataFrame
        Element-normalized populations for explicit levels, indexed by
        ``(atomic_number, ion_number, level_number)``.
    normalized_ion_populations : pandas.DataFrame
        Element-normalized total population for every retained ion stage,
        including the terminal bare-nucleus state.
    level_populations : pandas.DataFrame
        Absolute explicit-level populations reconstructed from the elemental
        number density.
    ion_populations : pandas.DataFrame
        Absolute total-ion populations reconstructed from the elemental
        number density.
    state_index : ElementalStateIndex
        Mapping between physical level/ion labels and matrix positions.
    """

    normalized_state_populations: pd.DataFrame
    normalized_level_populations: pd.DataFrame
    normalized_ion_populations: pd.DataFrame
    level_populations: pd.DataFrame
    ion_populations: pd.DataFrame
    state_index: ElementalStateIndex


class AnalyticEquilibriumIonPopulationSolver:
    """Solve ion populations using analytic radiative rates."""

    def __init__(
        self,
        rate_matrix_solver,
        photoionization_rate_solver,
        collisional_ionization_rate_solver,
        elemental_number_density,
        max_solver_iterations=100,
    ):
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver : EquilibriumIonRateMatrix
        photoionization_rate_solver : RadiativeIonizationRateSolver
        collisional_ionization_rate_solver : CollisionalIonizationRateSolver
        max_solver_iterations : int, optional
            Maximum number of iterations for the ion population solver, by default 100.
        """
        self.rate_matrix_solver = rate_matrix_solver
        self.photoionization_rate_solver = photoionization_rate_solver
        self.collisional_ionization_rate_solver = (
            collisional_ionization_rate_solver
        )
        self.elemental_number_density = elemental_number_density
        self.max_solver_iterations = max_solver_iterations

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
    ) -> ElementalPopulationSolution:
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
        ElementalPopulationSolution
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
            population = np.linalg.solve(
                normalized_rate_matrix,
                normalization,
            )
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
        return ElementalPopulationSolution(
            normalized_state_populations=state_populations,
            normalized_level_populations=normalized_level_populations,
            normalized_ion_populations=normalized_ion_populations,
            level_populations=level_populations,
            ion_populations=ion_populations,
            state_index=state_index,
        )

    def _solve_elemental_populations(
        self,
        rate_frames: tuple[pd.DataFrame, ...],
        fallback_ion_populations: pd.DataFrame,
    ) -> pd.DataFrame:
        """Solve all configured elements from extended rate matrices.

        Elements with rate data use the extended elemental matrix. Elements
        without configured continuum rates retain their current populations so
        their charge still contributes to the global closure.

        Parameters
        ----------
        rate_frames : tuple[pandas.DataFrame, ...]
            Photoionization, radiative recombination, collisional ionization,
            and collisional recombination rates.
        fallback_ion_populations : pandas.DataFrame
            Current absolute ion populations for elements without rate data.

        Returns
        -------
        pandas.DataFrame
            Absolute ion populations indexed by element and ion stage.
        """
        elemental_populations = []
        for atomic_number in self.elemental_number_density.index:
            elemental_rate_frames = tuple(
                rate_frame.loc[
                    rate_frame.index.get_level_values("atomic_number")
                    == atomic_number
                ]
                for rate_frame in rate_frames
            )
            if any(
                not rate_frame.empty for rate_frame in elemental_rate_frames
            ):
                matrix_set = self.rate_matrix_solver.solve_ion_and_level(
                    atomic_number,
                    photoion_rates_df=elemental_rate_frames[0],
                    recomb_rates_df=elemental_rate_frames[1],
                    collisional_ionization_rates_df=elemental_rate_frames[2],
                    collision_recombination_rates_df=elemental_rate_frames[3],
                )
                elemental_solution = self._build_elemental_population_solution(
                    matrix_set,
                    self.elemental_number_density.loc[atomic_number],
                )
                elemental_populations.append(elemental_solution.ion_populations)
            else:
                elemental_populations.append(
                    fallback_ion_populations.loc[[atomic_number]]
                )
        return pd.concat(elemental_populations)

    def _solve_trial_populations(
        self,
        electron_number_density: pd.Series,
        solve_trial_rates: Callable[[pd.Series], tuple[pd.DataFrame, ...]],
        fallback_ion_populations: pd.DataFrame,
    ) -> pd.DataFrame:
        """Solve elemental populations for one trial electron density.

        Parameters
        ----------
        electron_number_density : pandas.Series
            Trial electron densities indexed by shell.
        solve_trial_rates : callable
            Rate calculation for the supplied trial electron density.
        fallback_ion_populations : pandas.DataFrame
            Current ion populations for elements without rate data.

        Returns
        -------
        pandas.DataFrame
            Absolute ion populations for all elements.
        """
        return self._solve_elemental_populations(
            solve_trial_rates(electron_number_density),
            fallback_ion_populations,
        )

    def _solve_charge_conserving_populations(
        self,
        solve_trial_rates: Callable[[pd.Series], tuple[pd.DataFrame, ...]],
        fallback_ion_populations: pd.DataFrame,
    ) -> tuple[pd.Series, pd.DataFrame]:
        """Solve trial elemental populations and close their total charge."""
        solver = ChargeConservationSolver(
            self.elemental_number_density,
            partial(
                self._solve_trial_populations,
                solve_trial_rates=solve_trial_rates,
                fallback_ion_populations=fallback_ion_populations,
            ),
            max_solver_iterations=self.max_solver_iterations,
        )
        return solver.solve()

    def _calculate_trial_rates(
        self,
        trial_electron_number_density: pd.Series,
        base_electron_energy_distribution: ThermalElectronEnergyDistribution,
        radiation_field: object,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
    ) -> tuple[pd.DataFrame, ...]:
        """Calculate rates for a trial electron density.

        Parameters
        ----------
        trial_electron_number_density : pandas.Series
            Trial electron densities indexed by shell.
        base_electron_energy_distribution : ThermalElectronEnergyDistribution
            Energy and temperature used to construct the trial distribution.
        radiation_field : object
            Radiation field used by the photoionization solver.
        lte_level_population : pandas.DataFrame
            LTE level populations.
        estimated_level_population : pandas.DataFrame
            Estimated level populations.
        lte_ion_population : pandas.DataFrame
            LTE ion populations.
        estimated_ion_population : pandas.DataFrame
            Estimated ion populations.
        partition_function : pandas.DataFrame
            Partition function values.
        boltzmann_factor : pandas.DataFrame
            Boltzmann factors.

        Returns
        -------
        tuple[pandas.DataFrame, ...]
            Photoionization, recombination, collisional ionization, and
            collisional recombination rates.
        """
        trial_electron_distribution = ThermalElectronEnergyDistribution(
            base_electron_energy_distribution.energy,
            base_electron_energy_distribution.temperature,
            trial_electron_number_density.to_numpy() * u.cm**-3,
        )
        return self._calculate_rates(
            trial_electron_distribution,
            radiation_field,
            lte_level_population,
            estimated_level_population,
            lte_ion_population,
            estimated_ion_population,
            partition_function,
            boltzmann_factor,
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

    def solve(
        self,
        radiation_field: object,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        charge_conservation: bool = False,
        tolerance: float = 1e-14,
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Solves the normalized ion population values from the rate matrices.

        Parameters
        ----------
        radiation_field : RadiationField
            A radiation field that can compute its mean intensity.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        elemental_number_density : pd.DataFrame
            Elemental number density. Index is atomic number, columns are cells.
        lte_level_population : pd.DataFrame
            LTE level number density. Columns are cells.
        estimated_level_population : pd.DataFrame
            Estimated level number density. Columns are cells.
        lte_ion_population : pd.DataFrame
            LTE ion number density. Columns are cells.
        estimated_ion_population : pd.DataFrame
            Estimated ion number density. Columns are cells.
        charge_conservation : bool, optional
            Whether to include a charge conservation row in the rate matrix.
        tolerance : float, optional
            Tolerance for convergence of the ion population solver.

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
            if charge_conservation:
                calculate_trial_rates = partial(
                    self._calculate_trial_rates,
                    base_electron_energy_distribution=new_electron_energy_distribution,
                    radiation_field=radiation_field,
                    lte_level_population=lte_level_population,
                    estimated_level_population=estimated_level_population,
                    lte_ion_population=lte_ion_population,
                    estimated_ion_population=estimated_ion_population,
                    partition_function=partition_function,
                    boltzmann_factor=boltzmann_factor,
                )

                (
                    electron_population_solution,
                    ion_population_solution,
                ) = self._solve_charge_conserving_populations(
                    calculate_trial_rates,
                    estimated_ion_population,
                )
            else:
                (
                    photoion_rates_df,
                    recomb_rates_df,
                    collisional_ionization_rates_df,
                    collision_recombination_rates_df,
                ) = self._calculate_rates(
                    new_electron_energy_distribution,
                    radiation_field,
                    lte_level_population,
                    estimated_level_population,
                    lte_ion_population,
                    estimated_ion_population,
                    partition_function,
                    boltzmann_factor,
                )
                rates_matrices = self.rate_matrix_solver.solve(
                    photoion_rates_df,
                    recomb_rates_df,
                    collisional_ionization_rates_df,
                    collision_recombination_rates_df,
                    charge_conservation,
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
                        charge_conservation,
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
                    * ion_population_solution.index.get_level_values(
                        "ion_number"
                    )
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
                np.all(np.abs(delta_ion) < tolerance).any().any()
                and (np.abs(delta_electron) < tolerance).any().any()
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


class EstimatedEquilibriumIonPopulationSolver(
    AnalyticEquilibriumIonPopulationSolver
):
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


class LTEIonPopulationSolver(AnalyticEquilibriumIonPopulationSolver):
    """Solve ion populations from LTE rate matrices."""

    def __init__(
        self,
        rate_matrix_solver,
        elemental_number_density,
        max_solver_iterations=100,
    ):
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver : LTEIonRateMatrix
        elemental_number_density : pd.DataFrame
            Elemental number density. Index is atomic number, columns are cells.
        max_solver_iterations : int, optional
        """
        super().__init__(
            rate_matrix_solver,
            None,
            None,
            elemental_number_density,
            max_solver_iterations,
        )

    def solve(
        self,
        thermal_electron_energy_distribution,
        lte_ion_population,
        phi,
        partition_function,
        charge_conservation=False,
        tolerance=1e-14,
    ):
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
        charge_conservation : bool, optional
            Whether to include a charge conservation row in the rate matrix.
        tolerance : float, optional
            Tolerance for convergence of the ion population solver.

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
                charge_conservation,
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
                    charge_conservation,
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
                np.all(np.abs(delta_ion) < tolerance).any().any()
                and (np.abs(delta_electron) < tolerance).any().any()
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
