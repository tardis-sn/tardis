import logging

import astropy.units as u
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

LOWER_ION_LEVEL_H = 0


class IonPopulationSolver:
    def __init__(self, rate_matrix_solver, max_solver_iterations=100):
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver : IonRateMatrix
        """
        self.rate_matrix_solver = rate_matrix_solver
        self.max_solver_iterations = max_solver_iterations

    def __calculate_balance_vector(
        self,
        elemental_number_density,
        rate_matrix_index,
        charge_conservation=False,
    ):
        """Constructs the balance vector for the NLTE ionization solver set of
        equations by combining all solution vector blocks.

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
                elemental_number_density.loc[atomic_number],
            )
            for atomic_number in elemental_number_density.index
        ]

        if charge_conservation:
            balance_array.append(np.array([0.0]))

        return np.hstack(balance_array)

    def solve(
        self,
        radiation_field,
        thermal_electron_energy_distribution,
        elemental_number_density,
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
            number and ion number. Columns are cells.
        pd.DataFrame
            Normalized electron fraction values. Columns are cells.
        """
        # TODO: make more general indices that work for non-Hydrogen species
        # this is the i level in Lucy 2003
        lower_ion_level_index = (
            lte_level_population.index.get_level_values("ion_number")
            == LOWER_ION_LEVEL_H
        )

        # this is the k level in Lucy 2003
        upper_ion_population_index = (
            lte_ion_population.index.get_level_values("ion_number")
            > LOWER_ION_LEVEL_H
        )

        new_electron_energy_distribution = thermal_electron_energy_distribution

        for iteration in range(self.max_solver_iterations):
            self.rates_matrices = self.rate_matrix_solver.solve(
                radiation_field,
                new_electron_energy_distribution,
                lte_level_population.loc[lower_ion_level_index],
                estimated_level_population.loc[lower_ion_level_index],
                lte_ion_population.loc[upper_ion_population_index],
                estimated_ion_population.loc[upper_ion_population_index],
                partition_function,
                boltzmann_factor,
                charge_conservation,
            )
            solved_matrices = pd.DataFrame(
                index=self.rates_matrices.index,
                columns=self.rates_matrices.columns,
            )
            for cell in elemental_number_density.columns:
                balance_vector = self.__calculate_balance_vector(
                    elemental_number_density[cell],
                    self.rates_matrices.index,
                    charge_conservation,
                )
                solved_matrix = self.rates_matrices[cell].apply(
                    np.linalg.solve, args=(balance_vector,)
                )
                solved_matrices[cell] = solved_matrix

            ion_population_solution = pd.DataFrame(
                np.vstack(solved_matrices.values[0]).T,
                index=estimated_ion_population.index,
                columns=self.rates_matrices.columns,
            )

            if (ion_population_solution < 0).any().any():
                ion_population_solution[ion_population_solution < 0] = 0.0

            electron_population_solution = (
                ion_population_solution
                * ion_population_solution.index.get_level_values("ion_number")
                .values[np.newaxis]
                .T
            ).sum()

            delta_ion = (
                estimated_ion_population - ion_population_solution
            ) / ion_population_solution
            delta_electron = (
                new_electron_energy_distribution.number_density.value
                - electron_population_solution
            ) / electron_population_solution

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
