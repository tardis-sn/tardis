import logging

import astropy.units as u
import numpy as np
import numpy.typing as npt
import pandas as pd

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
    return level_population_at_lower_ion / (
        ion_population_at_upper_ion.values * electron_number_density.value
    )


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

    def solve(
        self,
        radiation_field,
        thermal_electron_energy_distribution,
        lte_level_population,
        estimated_level_population,
        lte_ion_population,
        estimated_ion_population,
        partition_function,
        boltzmann_factor,
        charge_conservation=False,
        tolerance=1e-14,
    ):
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
        # TODO: make more general indices that work for non-Hydrogen species
        # this is the i level in Lucy 2003
        lower_ion_level_index = get_lower_ion_level_index(lte_level_population)

        # this is the k level in Lucy 2003
        upper_ion_population_index = get_upper_ion_population_index(
            lte_ion_population
        )

        new_electron_energy_distribution = thermal_electron_energy_distribution

        for iteration in range(self.max_solver_iterations):
            photoion_rates_df, recomb_rates_df = (
                self.photoionization_rate_solver.solve(
                    radiation_field,
                    new_electron_energy_distribution,
                    lte_level_population.loc[lower_ion_level_index],
                    estimated_level_population.loc[lower_ion_level_index],
                    lte_ion_population.loc[upper_ion_population_index],
                    estimated_ion_population.loc[upper_ion_population_index],
                    partition_function,
                    boltzmann_factor,
                )
            )

            # Lucy 2003 Eq 14
            level_to_ion_population_factor = (
                calculate_level_to_ion_population_factor(
                    lte_level_population.loc[lower_ion_level_index],
                    lte_ion_population.loc[upper_ion_population_index],
                    new_electron_energy_distribution.number_density,
                )
            )

            (
                collisional_ionization_rates_df,
                collision_recombination_rates_df,
            ) = self.collisional_ionization_rate_solver.solve(
                new_electron_energy_distribution,
                level_to_ion_population_factor,
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

    def solve(
        self,
        estimated_radiation_field,
        thermal_electron_energy_distribution,
        lte_level_population,
        estimated_level_population,
        lte_ion_population,
        estimated_ion_population,
        partition_function,
        boltzmann_factor,
        charge_conservation=False,
        tolerance=1e-14,
    ):
        """Solves the normalized ion population values from the rate matrices.

        Parameters
        ----------
        estimated_radiation_field : dict
            Continuum radiation field estimators from Monte Carlo.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron properties.
        elemental_number_density : pd.DataFrame
            Elemental number density. Index is atomic number, columns are cells.
        time_simulation : u.Quantity
            Time of simulation.
        volume : u.Quantity
            Volume per cell.
        lte_level_population : pd.DataFrame
            LTE level number density. Columns are cells.
        estimated_level_population : pd.DataFrame
            Estimated level number density. Columns are cells.
        lte_ion_population : pd.DataFrame
            LTE ion number density. Columns are cells.
        estimated_ion_population : pd.DataFrame
            Estimated ion number density. Columns are cells.
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
        # TODO: make more general indices that work for non-Hydrogen species
        # this is the i level in Lucy 2003
        lower_ion_level_index = get_lower_ion_level_index(lte_level_population)

        # this is the k level in Lucy 2003
        upper_ion_population_index = get_upper_ion_population_index(
            lte_ion_population
        )

        new_electron_energy_distribution = thermal_electron_energy_distribution

        for iteration in range(self.max_solver_iterations):
            photoion_rates_df, recomb_rates_df = (
                self.photoionization_rate_solver.solve(
                    new_electron_energy_distribution,
                    estimated_radiation_field,
                    estimated_level_population.loc[lower_ion_level_index],
                    lte_level_population.loc[lower_ion_level_index],
                    estimated_ion_population.loc[upper_ion_population_index],
                    lte_ion_population.loc[upper_ion_population_index],
                )
            )

            # Lucy 2003 Eq 14
            level_to_ion_population_factor = (
                calculate_level_to_ion_population_factor(
                    lte_level_population.loc[lower_ion_level_index],
                    lte_ion_population.loc[upper_ion_population_index],
                    new_electron_energy_distribution.number_density,
                )
            )

            (
                collisional_ionization_rates_df,
                collision_recombination_rates_df,
            ) = self.collisional_ionization_rate_solver.solve(
                new_electron_energy_distribution,
                level_to_ion_population_factor,
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
