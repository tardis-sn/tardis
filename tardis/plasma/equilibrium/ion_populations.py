import logging

import astropy.units as u
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


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
        number_density : pandas.DataFrame
            Number densities of all present species.
        rate_matrix_index : pandas.MultiIndex
            (atomic_number, ion_number, treatment type)

        Returns
        -------
        numpy.array
        Solution vector for the NLTE ionization solver.
        """
        atomic_numbers = elemental_number_density.index
        balance_array = []
        for atomic_number in atomic_numbers:
            needed_number_of_elements = (
                rate_matrix_index.get_level_values("atomic_number")
                == atomic_number
            ).sum()
            balance_vector_block = np.zeros(needed_number_of_elements + 1)
            balance_vector_block[-1] = elemental_number_density.loc[
                atomic_number
            ]
            balance_array.append(balance_vector_block)
            if charge_conservation:
                balance_array.append(np.array([0.0]))
        return np.hstack(balance_array)

    def solve(
        self,
        radiation_field,
        thermal_electron_energy_distribution,
        elemental_number_density,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
        charge_conservation=False,
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

        new_electron_energy_distribution = thermal_electron_energy_distribution

        ion_population_solution = pd.DataFrame(
            index=ion_population.index, columns=ion_population.columns
        )

        for iteration in range(self.max_solver_iterations):
            self.rates_matrices = self.rate_matrix_solver.solve(
                radiation_field,
                new_electron_energy_distribution,
                lte_level_population.loc[lower_ion_level_index],
                level_population.loc[lower_ion_level_index],
                lte_ion_population.loc[upper_ion_population_index],
                ion_population.loc[upper_ion_population_index],
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
                index=ion_population.index,
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
                ion_population - ion_population_solution
            ) / ion_population_solution
            delta_electron = (
                new_electron_energy_distribution.number_density.value
                - electron_population_solution
            ) / electron_population_solution

            if (
                np.all(np.abs(delta_ion) < 1e-14).any().any()
                and (np.abs(delta_electron) < 1e-14).any().any()
            ):
                logger.info(
                    "Ion population solver converged after %d iterations.",
                    iteration + 1,
                )
                break

            ion_population = ion_population_solution
            new_electron_energy_distribution.number_density = (
                electron_population_solution.values * u.cm**-3
            )
        else:
            logger.warning(
                "Ion population solver did not converge after %d iterations.",
                iteration,
            )

        return ion_population_solution, electron_population_solution
