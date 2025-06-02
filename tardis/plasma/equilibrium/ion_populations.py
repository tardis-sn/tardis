import numpy as np
import pandas as pd


class IonPopulationSolver:
    max_solver_iterations = 100

    def __init__(self, rate_matrix_solver):
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver :
        """
        self.rate_matrix_solver = rate_matrix_solver

    def __calculate_ion_population(self, rates_matrix: np.ndarray):
        """Helper function to calculate the normalized, per-ion boltzmann factor.

        Parameters
        ----------
        rates_matrix : np.ndarray
            The rate matrix for a given species and cell.

        Returns
        -------
        np.ndarray
            The normalized, per-ion population.
        """
        normalized_ion_population = np.zeros(rates_matrix.shape[0])
        normalized_ion_population[1] = 1.0
        normalized_ion_population = np.linalg.solve(
            rates_matrix, normalized_ion_population
        )
        return normalized_ion_population[:-1]

    def solve(
        self,
        radiation_field,
        thermal_electron_energy_distribution,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        charge_conservation,
    ):
        """Solves the normalized ion population values from the rate matrices.

        Returns
        -------
        pd.DataFrame
            Normalized ion population values indexed by atomic number, ion
            number and ion number. Columns are cells.
        pd.DataFrame
            Normalized electron fraction values. Columns are cells.
        """
        # this is the i level in Lucy 2003
        lower_ion_level_index = (
            lte_level_population.index.get_level_values("ion_number") == 0
        )

        # this is the k level in Lucy 2003
        upper_ion_population_index = (
            lte_ion_population.index.get_level_values("ion_number") >= 1
        )

        self.rates_matrices = self.rate_matrix_solver.solve(
            radiation_field,
            thermal_electron_energy_distribution,
            lte_level_population.loc[lower_ion_level_index],
            level_population.loc[lower_ion_level_index],
            lte_ion_population.loc[upper_ion_population_index],
            ion_population.loc[upper_ion_population_index],
            charge_conservation,
        )

        iteration = 0
        electron_densities = 1

        while iteration < self.max_solver_iterations:
            normalized_ion_population = ion_population / ion_population.sum()

            self.rates_matrices = self.rate_matrix_solver.solve(
                radiation_field,
                thermal_electron_energy_distribution,
                lte_level_population.loc[lower_ion_level_index],
                level_population.loc[lower_ion_level_index],
                lte_ion_population.loc[upper_ion_population_index],
                ion_population.loc[upper_ion_population_index],
                charge_conservation,
            )
            solved_matrices = self.rates_matrices.map(
                self.__calculate_ion_population
            )

            ion_population_solution = pd.DataFrame(
                np.vstack(solved_matrices.values[0]).T,
                index=ion_population.index,
            )

            if (ion_population_solution < 0).any().any():
                ion_population_solution[ion_population_solution < 0] = 0.0

            electron_population_solution = ion_population_solution.sum()

            delta_ion = (
                normalized_ion_population - ion_population_solution
            ) / ion_population_solution
            delta_electron = (
                electron_densities - electron_population_solution
            ) / electron_population_solution

            if (
                np.all(np.abs(delta_ion) < 1e-8).any().any()
                and (np.abs(delta_electron) < 1e-8).any().any()
            ):
                break

            ion_population = ion_population_solution
            electron_densities = electron_population_solution

            iteration += 1

        return ion_population_solution, electron_population_solution
