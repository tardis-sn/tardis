import numpy as np
import pandas as pd


class IonPopulationSolver:
    max_solver_iterations = 100

    def __init__(self, rate_matrix_solver):
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver : IonRateMatrix
        """
        self.rate_matrix_solver = rate_matrix_solver

    def __calculate_ion_population(self, rates_matrix: np.ndarray):
        """Helper function to calculate the normalized ion population.

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
        # Number conservation, i.e. the sum of all ion populations is 1
        normalized_ion_population[1] = 1.0
        normalized_ion_population = np.linalg.solve(
            rates_matrix, normalized_ion_population
        )
        # drop charge conservation entry
        return normalized_ion_population[:-1]

    def solve(
        self,
        radiation_field,
        thermal_electron_energy_distribution,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        charge_conservation=True,
    ):
        """Solves the normalized ion population values from the rate matrices.

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
            Normalized ion population values indexed by atomic number, ion
            number and ion number. Columns are cells.
        pd.DataFrame
            Normalized electron fraction values. Columns are cells.
        """
        # TODO: make more general indices that work for non-Hydrogen species
        # this is the i level in Lucy 2003
        lower_ion_level_index = (
            lte_level_population.index.get_level_values("ion_number") == 0
        )

        # this is the k level in Lucy 2003
        upper_ion_population_index = (
            lte_ion_population.index.get_level_values("ion_number") >= 1
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

        print(iteration)

        return ion_population_solution, electron_population_solution
