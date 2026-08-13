import logging

import astropy.units as u
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.optimize import brentq

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix
from tardis.plasma.exceptions import PlasmaIonizationError
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
    PlanckianRadiationField,
)

logger = logging.getLogger(__name__)

LOWER_ION_LEVEL_H = 0
MINIMUM_BRENT_RELATIVE_TOLERANCE = 4 * np.finfo(float).eps
CHARGE_TOLERANCE = 1e-10


class IonPopulationSolver:
    """Solve ion populations from elemental ionization rate matrices."""

    def __init__(
        self,
        rate_matrix_solver: IonRateMatrix,
        max_solver_iterations: int = 100,
    ) -> None:
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver : IonRateMatrix
            Solver that builds ionization rate matrices.
        max_solver_iterations : int, optional
            Maximum iterations for lagged population-dependent corrections.
        """
        self.rate_matrix_solver = rate_matrix_solver
        self.max_solver_iterations = max_solver_iterations

    def _calculate_balance_vector(
        self,
        elemental_number_density: pd.Series,
        rate_matrix_index: pd.MultiIndex,
        charge_conservation: bool = False,
    ) -> npt.NDArray[np.float64]:
        """Construct the balance vector for the NLTE ionization solver.

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

    def solve_element_populations_at_electron_density(
        self,
        electron_density: npt.NDArray[np.float64],
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        elemental_number_density: pd.DataFrame,
    ) -> pd.DataFrame:
        """Solve elemental ion populations at supplied electron densities.

        Parameters
        ----------
        electron_density : numpy.ndarray
            Electron number densities in cm^-3 ordered like the output columns.
        radiation_field : RadiationField
            Fixed radiation field.
        thermal_electron_energy_distribution : tardis.plasma.electron_energy_distribution.ThermalElectronEnergyDistribution
            Fixed electron energy and temperature inputs.
        lte_level_population : pandas.DataFrame
            LTE level number density. Columns are cells.
        estimated_level_population : pandas.DataFrame
            Lagged estimated level number density. Columns are cells.
        lte_ion_population : pandas.DataFrame
            LTE ion number density. Columns are cells.
        estimated_ion_population : pandas.DataFrame
            Lagged estimated ion number density. Columns are cells.
        elemental_number_density : pandas.DataFrame
            Elemental number density. Index is atomic number, columns are cells.

        Returns
        -------
        pandas.DataFrame
            Absolute ion populations indexed by atomic number and ion number.
        """
        trial_electron_distribution = ThermalElectronEnergyDistribution(
            thermal_electron_energy_distribution.energy,
            thermal_electron_energy_distribution.temperature,
            electron_density * u.cm**-3,
        )
        self.rates_matrices = self.rate_matrix_solver.solve(
            radiation_field,
            trial_electron_distribution,
            lte_level_population,
            estimated_level_population,
            lte_ion_population,
            estimated_ion_population,
            partition_function,
            boltzmann_factor,
        )
        ion_population_index = self.rate_matrix_solver.ion_population_index
        ion_population = pd.DataFrame(
            0.0,
            index=ion_population_index,
            columns=elemental_number_density.columns,
        )
        for shell in ion_population.columns:
            for atomic_number in self.rates_matrices.index:
                matrix = self.rates_matrices.loc[atomic_number, shell]
                if not np.all(np.isfinite(matrix)):
                    raise PlasmaIonizationError(
                        "Nonfinite ion population matrix for atomic number "
                        f"{atomic_number}, shell {shell}."
                    )
                right_hand_side = np.zeros(matrix.shape[0])
                right_hand_side[1] = 1.0
                try:
                    normalized_population = np.linalg.solve(
                        matrix, right_hand_side
                    )
                except np.linalg.LinAlgError as exc:
                    raise PlasmaIonizationError(
                        "Singular ion population matrix for atomic number "
                        f"{atomic_number}, shell {shell}."
                    ) from exc
                if not np.all(np.isfinite(normalized_population)):
                    raise PlasmaIonizationError(
                        "Nonfinite ion population for atomic number "
                        f"{atomic_number}, shell {shell}."
                    )
                minimum_population = normalized_population.min()
                if minimum_population < -1e-12:
                    raise PlasmaIonizationError(
                        "Negative ion population for atomic number "
                        f"{atomic_number}, shell {shell}: "
                        f"{minimum_population}."
                    )
                normalized_population[normalized_population < 0.0] = 0.0
                normalized_population /= normalized_population.sum()
                population_index = ion_population_index[
                    ion_population_index.get_level_values("atomic_number")
                    == atomic_number
                ]
                if len(population_index) != matrix.shape[0]:
                    raise PlasmaIonizationError(
                        "Ion population index size does not match matrix "
                        f"shape for atomic number {atomic_number}, "
                        f"shell {shell}."
                    )
                ion_population.loc[population_index, shell] = (
                    normalized_population
                    * elemental_number_density.loc[atomic_number, shell]
                )
        return ion_population

    def solve_charge_balance(
        self,
        electron_density: npt.NDArray[np.float64],
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        elemental_number_density: pd.DataFrame,
        maximum_electron_density: pd.Series,
    ) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
        """Solve all shells and return their normalized charge residuals."""
        ion_population = self.solve_element_populations_at_electron_density(
            electron_density,
            radiation_field,
            thermal_electron_energy_distribution,
            lte_level_population,
            estimated_level_population,
            lte_ion_population,
            estimated_ion_population,
            partition_function,
            boltzmann_factor,
            elemental_number_density,
        )
        electron_population = pd.Series(
            electron_density, index=elemental_number_density.columns
        )
        charge_density = (
            ion_population
            * ion_population.index.get_level_values("ion_number").to_numpy()[:, None]
        ).sum()
        charge_residual = (charge_density - electron_population) / (
            maximum_electron_density.replace(0.0, 1.0)
        )
        return ion_population, electron_population, charge_residual

    def solve_shell_charge(
        self,
        shell_idx: int,
        maximum_electron_density: float,
        base_electron_density: npt.NDArray[np.float64],
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        elemental_number_density: pd.DataFrame,
        maximum_electron_densities: pd.Series,
    ) -> float:
        """Solve the scalar charge equation for one shell."""
        if maximum_electron_density == 0.0:
            return 0.0

        def charge_residual(
            electron_density_fraction: float,
        ) -> float:
            electron_density = base_electron_density.copy()
            electron_density[shell_idx] = (
                electron_density_fraction * maximum_electron_density
            )
            return self.solve_charge_balance(
                electron_density,
                radiation_field,
                thermal_electron_energy_distribution,
                lte_level_population,
                estimated_level_population,
                lte_ion_population,
                estimated_ion_population,
                partition_function,
                boltzmann_factor,
                elemental_number_density,
                maximum_electron_densities,
            )[2].iloc[shell_idx]

        lower_residual = charge_residual(0.0)
        upper_residual = charge_residual(1.0)
        if lower_residual == 0.0:
            return 0.0
        if upper_residual == 0.0:
            return maximum_electron_density
        if lower_residual * upper_residual > 0.0:
            raise PlasmaIonizationError(
                f"Charge residual does not bracket shell {shell_idx}: "
                f"Q_hat(0)={lower_residual}, Q_hat(1)={upper_residual}."
            )
        electron_density_fraction = brentq(
            charge_residual,
            0.0,
            1.0,
            xtol=CHARGE_TOLERANCE,
            rtol=MINIMUM_BRENT_RELATIVE_TOLERANCE,
            maxiter=self.max_solver_iterations,
        )
        return electron_density_fraction * maximum_electron_density

    def solve_charge_conserving(
        self,
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        elemental_number_density: pd.DataFrame,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame | float,
        boltzmann_factor: pd.DataFrame,
        tolerance: float,
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Solve ion populations with one shared electron density per shell."""
        electron_density = (
            thermal_electron_energy_distribution.number_density.cgs.value.copy()
        )
        atomic_numbers = elemental_number_density.index.to_numpy(dtype=float)
        maximum_electron_density = elemental_number_density.multiply(
            atomic_numbers, axis=0
        ).sum()
        for iteration in range(self.max_solver_iterations):
            for shell_idx in range(len(elemental_number_density.columns)):
                electron_density[shell_idx] = self.solve_shell_charge(
                    shell_idx,
                    maximum_electron_density.iloc[shell_idx],
                    electron_density,
                    radiation_field,
                    thermal_electron_energy_distribution,
                    lte_level_population,
                    estimated_level_population,
                    lte_ion_population,
                    estimated_ion_population,
                    partition_function,
                    boltzmann_factor,
                    elemental_number_density,
                    maximum_electron_density,
                )
            (
                ion_population_solution,
                electron_population_solution,
                charge_residual,
            ) = self.solve_charge_balance(
                electron_density,
                radiation_field,
                thermal_electron_energy_distribution,
                lte_level_population,
                estimated_level_population,
                lte_ion_population,
                estimated_ion_population,
                partition_function,
                boltzmann_factor,
                elemental_number_density,
                maximum_electron_density,
            )
            common_index = estimated_ion_population.index.intersection(
                ion_population_solution.index
            )
            if len(common_index) == 0:
                population_converged = True
            else:
                population_delta = (
                    estimated_ion_population.loc[common_index]
                    - ion_population_solution.loc[common_index]
                ) / np.maximum(
                    np.abs(ion_population_solution.loc[common_index]), 1e-300
                )
                population_converged = np.all(
                    np.abs(population_delta) < tolerance
                ).all()
            if (
                np.all(np.abs(charge_residual) < CHARGE_TOLERANCE)
                and population_converged
            ):
                logger.info(
                    "Ion population solver converged after %d iterations.",
                    iteration + 1,
                )
                return ion_population_solution, electron_population_solution
            if len(common_index) != 0:
                estimated_ion_population = ion_population_solution.loc[
                    estimated_ion_population.index
                ]

        raise PlasmaIonizationError(
            "Ion population solver did not converge after "
            f"{self.max_solver_iterations} iterations."
        )

    def solve_non_charge_conserving(
        self,
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        elemental_number_density: pd.DataFrame,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        tolerance: float,
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Solve ion populations without imposing charge conservation.

        Parameters
        ----------
        radiation_field : DilutePlanckianRadiationField or PlanckianRadiationField
            Radiation field used to calculate photoionization rates.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution, including the current electron
            number density.
        elemental_number_density : pandas.DataFrame
            Elemental number density indexed by atomic number, with cells in
            the columns.
        lte_level_population : pandas.DataFrame
            LTE level number density.
        estimated_level_population : pandas.DataFrame
            Current estimated level number density.
        lte_ion_population : pandas.DataFrame
            LTE ion number density.
        estimated_ion_population : pandas.DataFrame
            Current estimated ion number density.
        partition_function : pandas.DataFrame or float
            Partition function used by the rate solvers.
        boltzmann_factor : pandas.DataFrame
            Level Boltzmann factors.
        tolerance : float
            Relative convergence tolerance for ion and electron populations.

        Returns
        -------
        tuple of pandas.DataFrame and pandas.Series
            Ion population and electron number density solutions.
        """
        new_electron_energy_distribution = thermal_electron_energy_distribution

        for iteration in range(self.max_solver_iterations):
            self.rates_matrices = self.rate_matrix_solver.solve(
                radiation_field,
                new_electron_energy_distribution,
                lte_level_population,
                estimated_level_population,
                lte_ion_population,
                estimated_ion_population,
                partition_function,
                boltzmann_factor,
            )
            solved_matrices = pd.DataFrame(
                index=self.rates_matrices.index,
                columns=self.rates_matrices.columns,
            )
            for cell in elemental_number_density.columns:
                balance_vector = self._calculate_balance_vector(
                    elemental_number_density[cell],
                    self.rates_matrices.index,
                )
                solved_matrix = self.rates_matrices[cell].apply(
                    np.linalg.solve, args=(balance_vector,)
                )
                solved_matrices[cell] = solved_matrix

            ion_population_solution = pd.DataFrame(
                np.vstack(solved_matrices.values[0]).T,
                index=self.rate_matrix_solver.ion_population_index,
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

            estimated_solution = ion_population_solution.loc[
                estimated_ion_population.index
            ]
            delta_ion = (
                estimated_ion_population - estimated_solution
            ) / estimated_solution
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

            estimated_ion_population = estimated_solution
            new_electron_energy_distribution.number_density = (
                electron_population_solution.values * u.cm**-3
            )
        else:
            logger.warning(
                "Ion population solver did not converge after %d iterations.",
                iteration,
            )

        return ion_population_solution, electron_population_solution

    def solve(
        self,
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        elemental_number_density: pd.DataFrame,
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
        partition_function : pandas.DataFrame or float
            Partition function values used by the rate solver.
        boltzmann_factor : pd.DataFrame
            Boltzmann factors used by the rate solver.
        charge_conservation : bool, optional
            Whether to solve one shared charge-conserving electron density per
            shell.
        tolerance : float, optional
            Tolerance for convergence of the ion population solver.

        Returns
        -------
        pd.DataFrame
            Ion population values indexed by atomic number and ion number.
            Columns are cells.
        pd.Series
            Electron number densities indexed by cell.
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

        if charge_conservation:
            return self.solve_charge_conserving(
                radiation_field,
                thermal_electron_energy_distribution,
                elemental_number_density,
                lte_level_population.loc[lower_ion_level_index],
                estimated_level_population.loc[lower_ion_level_index],
                lte_ion_population.loc[upper_ion_population_index],
                estimated_ion_population.loc[upper_ion_population_index],
                partition_function,
                boltzmann_factor,
                tolerance,
            )

        return self.solve_non_charge_conserving(
            radiation_field,
            thermal_electron_energy_distribution,
            elemental_number_density,
            lte_level_population.loc[lower_ion_level_index],
            estimated_level_population.loc[lower_ion_level_index],
            lte_ion_population.loc[upper_ion_population_index],
            estimated_ion_population.loc[upper_ion_population_index],
            partition_function,
            boltzmann_factor,
            tolerance,
        )
