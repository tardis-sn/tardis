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

MINIMUM_BRENT_RELATIVE_TOLERANCE = 4 * np.finfo(float).eps
CHARGE_TOLERANCE = 1e-10
MINIMUM_ELECTRON_DENSITY_FRACTION = np.finfo(float).eps


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
        level_to_continuum_saha_factor: pd.DataFrame,
        elemental_number_density: pd.DataFrame,
    ) -> npt.NDArray[np.float64]:
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
            Previous estimated level number density. Columns are cells.
        lte_ion_population : pandas.DataFrame
            LTE ion number density. Columns are cells.
        estimated_ion_population : pandas.DataFrame
            Previous estimated ion number density. Columns are cells.
        level_to_continuum_saha_factor : pandas.DataFrame
            Density-independent Lucy level-to-continuum Saha factor.
        elemental_number_density : pandas.DataFrame
            Elemental number density. Index is atomic number, columns are cells.

        Returns
        -------
        numpy.ndarray
            Absolute ion populations ordered like ``ion_population_index`` and
            the elemental-density columns.
        """
        rate_matrices = self.rate_matrix_solver.solve(
            radiation_field,
            thermal_electron_energy_distribution,
            lte_level_population,
            estimated_level_population,
            lte_ion_population,
            estimated_ion_population,
            partition_function,
            boltzmann_factor,
            level_to_continuum_saha_factor,
        )

        ion_population_index = self.rate_matrix_solver.ion_population_index
        ion_population = np.zeros(
            (len(ion_population_index), len(rate_matrices.columns))
        )

        atomic_numbers = ion_population_index.get_level_values(
            "atomic_number"
        ).to_numpy()

        rate_matrix_atomic_numbers = rate_matrices.index.to_numpy()
        rate_matrix_arrays = rate_matrices.to_numpy()

        population_indices = [
            np.flatnonzero(atomic_numbers == atomic_number)
            for atomic_number in rate_matrix_atomic_numbers
        ]

        elemental_indices = elemental_number_density.index.get_indexer(
            rate_matrix_atomic_numbers
        )
        elemental_number_density_array = elemental_number_density.to_numpy()

        for atomic_number_idx, atomic_number in enumerate(
            rate_matrix_atomic_numbers
        ):
            matrices = np.stack(rate_matrix_arrays[atomic_number_idx])
            nonfinite_matrices = ~np.isfinite(matrices).all(axis=(1, 2))

            if np.any(nonfinite_matrices):
                shell_idx = np.flatnonzero(nonfinite_matrices)[0]
                raise PlasmaIonizationError(
                    "Nonfinite ion population matrix for atomic number "
                    f"{atomic_number}, shell {rate_matrices.columns[shell_idx]}."
                )

            right_hand_side = np.zeros((len(matrices), matrices.shape[1], 1))
            right_hand_side[:, 1] = 1.0

            normalized_population = np.linalg.solve(matrices, right_hand_side)[
                :, :, 0
            ]

            nonfinite_populations = ~np.isfinite(normalized_population).all(
                axis=1
            )
            minimum_population = normalized_population.min(axis=1)
            if np.any(nonfinite_populations) or np.any(
                minimum_population < -1e-12
            ):
                shell_idx = np.flatnonzero(nonfinite_populations)[0]
                raise PlasmaIonizationError(
                    "Nonfinite or negative ion population for atomic number "
                    f"{atomic_number}, shell {rate_matrices.columns[shell_idx]}."
                    f"{minimum_population[shell_idx]}."
                )

            normalized_population[normalized_population < 0.0] = 0.0
            normalized_population /= normalized_population.sum(axis=1)[:, None]
            population_idx = population_indices[atomic_number_idx]

            elemental_idx = elemental_indices[atomic_number_idx]
            ion_population[population_idx] = (
                normalized_population.T
                * elemental_number_density_array[elemental_idx]
            )

        self.rates_matrices = rate_matrices
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
        level_to_continuum_saha_factor: pd.DataFrame,
        elemental_number_density: pd.DataFrame,
        maximum_electron_density: npt.NDArray[np.float64],
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Solve ion populations and calculate the normalized charge residual.

        Parameters
        ----------
        electron_density : npt.NDArray[np.float64]
            Trial electron number densities for each shell.
        radiation_field : DilutePlanckianRadiationField | PlanckianRadiationField
            Radiation field used to calculate ionization rates.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution used by the rate-matrix solver.
        lte_level_population : pd.DataFrame
            LTE level populations used by the rate-matrix solver.
        estimated_level_population : pd.DataFrame
            Estimated level populations used for lagged rate corrections.
        lte_ion_population : pd.DataFrame
            LTE ion populations used by the rate-matrix solver.
        estimated_ion_population : pd.DataFrame
            Estimated ion populations used for lagged rate corrections.
        partition_function : pd.DataFrame
            Partition functions used by the rate-matrix solver.
        boltzmann_factor : pd.DataFrame
            Boltzmann factors used by the rate-matrix solver.
        level_to_continuum_saha_factor : pd.DataFrame
            Density-independent Lucy level-to-continuum Saha factor.
        elemental_number_density : pd.DataFrame
            Elemental number densities indexed by atomic number and shell.
        maximum_electron_density : npt.NDArray[np.float64]
            Maximum possible electron number density for each shell, used to
            normalize the charge residual.

        Returns
        -------
        tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]
            Absolute ion populations and normalized charge residuals for each
            shell.
        """
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
            level_to_continuum_saha_factor,
            elemental_number_density,
        )

        charge_density = (
            ion_population
            * self.rate_matrix_solver.ion_population_index.get_level_values(
                "ion_number"
            ).to_numpy()[:, None]
        ).sum(axis=0)

        charge_residual = (charge_density - electron_density) / np.where(
            maximum_electron_density == 0.0, 1.0, maximum_electron_density
        )
        return ion_population, charge_residual

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
        level_to_continuum_saha_factor: pd.DataFrame,
        elemental_number_density: pd.DataFrame,
        maximum_electron_densities: npt.NDArray[np.float64],
    ) -> float:
        """Solve the charge balance for one shell.

        Parameters
        ----------
        shell_idx : int
            Index of the shell whose electron density is being solved.
        maximum_electron_density : float
            Maximum possible electron number density in the shell.
        base_electron_density : npt.NDArray[np.float64]
            Electron number densities used for all shells before updating the
            selected shell.
        radiation_field : DilutePlanckianRadiationField | PlanckianRadiationField
            Radiation field used to calculate ionization rates.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution used by the rate-matrix solver.
        lte_level_population : pd.DataFrame
            LTE level populations used by the rate-matrix solver.
        estimated_level_population : pd.DataFrame
            Estimated level populations used for lagged rate corrections.
        lte_ion_population : pd.DataFrame
            LTE ion populations used by the rate-matrix solver.
        estimated_ion_population : pd.DataFrame
            Estimated ion populations used for lagged rate corrections.
        partition_function : pd.DataFrame
            Partition functions used by the rate-matrix solver.
        boltzmann_factor : pd.DataFrame
            Boltzmann factors used by the rate-matrix solver.
        level_to_continuum_saha_factor : pd.DataFrame
            Density-independent Lucy level-to-continuum Saha factor.
        elemental_number_density : pd.DataFrame
            Elemental number densities indexed by atomic number and shell.
        maximum_electron_densities : npt.NDArray[np.float64]
            Maximum possible electron number density for each shell, used to
            normalize the charge residual.

        Returns
        -------
        float
            Charge-balanced electron number density for the selected shell.

        Raises
        ------
        PlasmaIonizationError
            If the charge residual is not bracketed over the allowed electron
            density interval.
        """

        def charge_residual(
            electron_density_fraction: float,
        ) -> float:
            """Calculate the normalized charge residual for one trial density."""
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
                level_to_continuum_saha_factor,
                elemental_number_density,
                maximum_electron_densities,
            )[1][shell_idx]

        try:
            electron_density_fraction = brentq(
                charge_residual,
                MINIMUM_ELECTRON_DENSITY_FRACTION,
                1.0,
                xtol=CHARGE_TOLERANCE,
                rtol=MINIMUM_BRENT_RELATIVE_TOLERANCE,
                maxiter=self.max_solver_iterations,
            )
        except ValueError as exc:
            lower_residual = charge_residual(MINIMUM_ELECTRON_DENSITY_FRACTION)
            upper_residual = charge_residual(1.0)
            raise PlasmaIonizationError(
                f"Charge residual does not bracket shell {shell_idx}: "
                "Q_hat(near-neutral)="
                f"{lower_residual}, Q_hat(1)={upper_residual}."
            ) from exc

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
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        level_to_continuum_saha_factor: pd.DataFrame,
        tolerance: float,
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Solve ion populations while enforcing charge conservation.

        Parameters
        ----------
        radiation_field : DilutePlanckianRadiationField | PlanckianRadiationField
            Radiation field used to calculate ionization rates.
        thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
            Electron energy distribution, including the initial electron
            number densities.
        elemental_number_density : pd.DataFrame
            Elemental number densities indexed by atomic number and shell.
        lte_level_population : pd.DataFrame
            LTE level populations used by the rate-matrix solver.
        estimated_level_population : pd.DataFrame
            Estimated level populations used for lagged rate corrections.
        lte_ion_population : pd.DataFrame
            LTE ion populations used by the rate-matrix solver.
        estimated_ion_population : pd.DataFrame
            Estimated ion populations used for lagged rate corrections.
        partition_function : pd.DataFrame
            Partition functions used by the rate-matrix solver.
        boltzmann_factor : pd.DataFrame
            Boltzmann factors used by the rate-matrix solver.
        level_to_continuum_saha_factor : pd.DataFrame
            Density-independent Lucy level-to-continuum Saha factor.
        tolerance : float
            Relative convergence tolerance for the ion populations.

        Returns
        -------
        tuple[pd.DataFrame, pd.Series]
            Ion populations indexed by atomic number and ion number, and
            charge-balanced electron number densities indexed by shell.

        Raises
        ------
        PlasmaIonizationError
            If the ion population solver does not converge within the maximum
            number of iterations or a shell charge balance cannot be bracketed.
        """
        electron_density = (
            thermal_electron_energy_distribution.number_density.cgs.value.copy()
        )

        atomic_numbers = elemental_number_density.index.to_numpy(dtype=float)

        maximum_electron_density = elemental_number_density.multiply(
            atomic_numbers, axis=0
        ).sum()
        maximum_electron_density_array = maximum_electron_density.to_numpy()

        converged_shells = np.zeros(
            len(elemental_number_density.columns), dtype=bool
        )

        estimated_population_indices: npt.NDArray[np.intp] | None = None
        solution_population_indices: npt.NDArray[np.intp] | None = None

        for iteration in range(self.max_solver_iterations):
            logger.info("Ion solver iteration %d", iteration + 1)
            for shell_idx in np.flatnonzero(~converged_shells):
                electron_density[shell_idx] = self.solve_shell_charge(
                    shell_idx,
                    maximum_electron_density_array[shell_idx],
                    electron_density,
                    radiation_field,
                    thermal_electron_energy_distribution,
                    lte_level_population,
                    estimated_level_population,
                    lte_ion_population,
                    estimated_ion_population,
                    partition_function,
                    boltzmann_factor,
                    level_to_continuum_saha_factor,
                    elemental_number_density,
                    maximum_electron_density_array,
                )
            ion_population_solution, charge_residual = (
                self.solve_charge_balance(
                    electron_density,
                    radiation_field,
                    thermal_electron_energy_distribution,
                    lte_level_population,
                    estimated_level_population,
                    lte_ion_population,
                    estimated_ion_population,
                    partition_function,
                    boltzmann_factor,
                    level_to_continuum_saha_factor,
                    elemental_number_density,
                    maximum_electron_density_array,
                )
            )
            if estimated_population_indices is None:
                solution_indices = (
                    self.rate_matrix_solver.ion_population_index.get_indexer(
                        estimated_ion_population.index
                    )
                )
                estimated_population_indices = np.flatnonzero(
                    solution_indices >= 0
                )
                solution_population_indices = solution_indices[
                    estimated_population_indices
                ]
            if len(estimated_population_indices) == 0:
                population_converged = np.ones_like(
                    converged_shells, dtype=bool
                )
            else:
                population_delta = (
                    estimated_ion_population.to_numpy()[
                        estimated_population_indices
                    ]
                    - ion_population_solution[solution_population_indices]
                ) / np.maximum(
                    np.abs(
                        ion_population_solution[solution_population_indices]
                    ),
                    1e-300,
                )
                population_converged = np.all(
                    np.abs(population_delta) < tolerance, axis=0
                )

            converged_shells |= (
                np.abs(charge_residual) < CHARGE_TOLERANCE
            ) & population_converged
            if np.all(converged_shells):
                logger.info(
                    "Ion population solver converged after %d iterations.",
                    iteration + 1,
                )
                return (
                    pd.DataFrame(
                        ion_population_solution,
                        index=self.rate_matrix_solver.ion_population_index,
                        columns=elemental_number_density.columns,
                    ),
                    pd.Series(
                        electron_density,
                        index=elemental_number_density.columns,
                    ),
                )

            if len(estimated_population_indices) != 0:
                estimated_population_values = estimated_ion_population.to_numpy(
                    copy=True
                )
                estimated_population_values[estimated_population_indices] = (
                    ion_population_solution[solution_population_indices]
                )
                estimated_ion_population = pd.DataFrame(
                    estimated_population_values,
                    index=estimated_ion_population.index,
                    columns=estimated_ion_population.columns,
                )

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
        level_to_continuum_saha_factor: pd.DataFrame | None = None,
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
        level_to_continuum_saha_factor : pd.DataFrame, optional
            Density-independent Lucy level-to-continuum Saha factor. Required
            when charge conservation varies the trial electron density.

        Returns
        -------
        pd.DataFrame
            Ion population values indexed by atomic number and ion number.
            Columns are cells.
        pd.Series
            Electron number densities indexed by cell.
        """
        bound_level_index = lte_level_population.index.get_level_values(
            "ion_number"
        ) < lte_level_population.index.get_level_values("atomic_number")
        continuum_ion_index = (
            lte_ion_population.index.get_level_values("ion_number") > 0
        )

        if charge_conservation:
            if level_to_continuum_saha_factor is None:
                raise ValueError(
                    "level_to_continuum_saha_factor is required when "
                    "charge_conservation is enabled."
                )
            return self.solve_charge_conserving(
                radiation_field,
                thermal_electron_energy_distribution,
                elemental_number_density,
                lte_level_population.loc[bound_level_index],
                estimated_level_population.loc[bound_level_index],
                lte_ion_population.loc[continuum_ion_index],
                estimated_ion_population.loc[continuum_ion_index],
                partition_function,
                boltzmann_factor.loc[bound_level_index],
                level_to_continuum_saha_factor.loc[
                    lte_level_population.index[bound_level_index]
                ],
                tolerance,
            )

        return self.solve_non_charge_conserving(
            radiation_field,
            thermal_electron_energy_distribution,
            elemental_number_density,
            lte_level_population.loc[bound_level_index],
            estimated_level_population.loc[bound_level_index],
            lte_ion_population.loc[continuum_ion_index],
            estimated_ion_population.loc[continuum_ion_index],
            partition_function,
            boltzmann_factor.loc[bound_level_index],
            tolerance,
        )
