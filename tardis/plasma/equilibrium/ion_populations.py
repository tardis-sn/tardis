import logging
from time import perf_counter

import astropy.units as u
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.optimize import brentq

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rate_matrix import (
    AnalyticIonRateMatrix,
    EstimatedIonRateMatrix,
)
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
        rate_matrix_solver: AnalyticIonRateMatrix | EstimatedIonRateMatrix,
        max_solver_iterations: int = 100,
    ) -> None:
        """Solve the normalized ion population values from the rate matrices.

        Parameters
        ----------
        rate_matrix_solver : AnalyticIonRateMatrix | EstimatedIonRateMatrix
            Solver that builds ionization rate matrices.
        max_solver_iterations : int, optional
            Maximum iterations for lagged population-dependent corrections.
        """
        self.rate_matrix_solver = rate_matrix_solver
        self.max_solver_iterations = max_solver_iterations
        self._logged_frozen_ion_populations: set[tuple[int, object]] = set()
        self._logged_negligible_ion_populations: set[tuple[int, object]] = set()
        self._charge_root_diagnostics: dict[
            object, tuple[float, int, int, bool, str]
        ] = {}

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
        | PlanckianRadiationField
        | None,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        level_to_continuum_saha_factor: pd.DataFrame,
        elemental_number_density: pd.DataFrame,
        lte_ionization_factor: pd.DataFrame | None = None,
        maximum_electron_density: npt.NDArray[np.float64] | None = None,
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
        partition_function : pandas.DataFrame
            Partition functions used by the rate-matrix solver.
        boltzmann_factor : pandas.DataFrame
            Boltzmann factors used by the rate-matrix solver.
        elemental_number_density : pandas.DataFrame
            Elemental number density. Index is atomic number, columns are cells.
        lte_ionization_factor : pandas.DataFrame, optional
            LTE ionization factors used to normalize ionization rates.
        maximum_electron_density : numpy.ndarray, optional
            Maximum possible electron density in each shell. When supplied by
            the charge-balance solve, it bounds the charge contribution of an
            element whose singular rate matrix has no usable prior population.

        Returns
        -------
        numpy.ndarray
            Absolute ion populations ordered like ``ion_population_index`` and
            the elemental-density columns.
        """
        if (
            getattr(
                self.rate_matrix_solver,
                "prepared_matrix_coefficients",
                None,
            )
            is None
        ):
            trial_electron_distribution = ThermalElectronEnergyDistribution(
                thermal_electron_energy_distribution.energy,
                thermal_electron_energy_distribution.temperature,
                electron_density * u.cm**-3,
            )
            rate_matrices = self.rate_matrix_solver.solve(
                radiation_field,
                trial_electron_distribution,
                lte_level_population,
                estimated_level_population,
                lte_ion_population,
                estimated_ion_population,
                partition_function,
                boltzmann_factor,
                level_to_continuum_saha_factor,
                lte_ionization_factor=lte_ionization_factor,
            )
        else:
            rate_matrices = self.rate_matrix_solver.solve_prepared(
                electron_density,
                lte_level_population.columns,
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
        aligned_estimated_ion_population = estimated_ion_population.reindex(
            index=ion_population_index,
            columns=rate_matrices.columns,
            fill_value=0.0,
        ).fillna(0.0)

        for atomic_number_idx, atomic_number in enumerate(
            rate_matrix_atomic_numbers
        ):
            elemental_idx = elemental_indices[atomic_number_idx]
            if elemental_idx < 0:
                continue
            matrices = np.stack(rate_matrix_arrays[atomic_number_idx])
            elemental_density = elemental_number_density_array[elemental_idx]
            populated_shell_indices = np.flatnonzero(elemental_density != 0.0)
            population_idx = population_indices[atomic_number_idx]
            if len(populated_shell_indices) == 0:
                continue

            populated_matrices = matrices[populated_shell_indices]
            nonfinite_matrices = ~np.isfinite(populated_matrices).all(
                axis=(1, 2)
            )

            if np.any(nonfinite_matrices):
                populated_shell_idx = np.flatnonzero(nonfinite_matrices)[0]
                shell_idx = populated_shell_indices[populated_shell_idx]
                raise PlasmaIonizationError(
                    "Nonfinite ion population matrix for atomic number "
                    f"{atomic_number}, shell {rate_matrices.columns[shell_idx]}."
                )

            right_hand_side = np.zeros(
                (
                    len(populated_matrices),
                    populated_matrices.shape[1],
                    1,
                )
            )
            right_hand_side[:, 1] = 1.0
            populated_normalized_population = np.zeros(
                (len(populated_matrices), populated_matrices.shape[1]),
                dtype=np.float64,
            )
            skipped_populated_shells = np.zeros(
                len(populated_matrices), dtype=bool
            )

            try:
                populated_normalized_population = np.linalg.solve(
                    populated_matrices, right_hand_side
                )[:, :, 0]
            except np.linalg.LinAlgError as exc:
                for populated_shell_idx, (
                    matrix,
                    shell_right_hand_side,
                ) in enumerate(
                    zip(populated_matrices, right_hand_side, strict=True)
                ):
                    try:
                        populated_normalized_population[populated_shell_idx] = (
                            np.linalg.solve(matrix, shell_right_hand_side)[:, 0]
                        )
                    except np.linalg.LinAlgError:
                        shell_idx = populated_shell_indices[populated_shell_idx]
                        shell = rate_matrices.columns[shell_idx]
                        estimated_population = (
                            aligned_estimated_ion_population.loc[
                                ion_population_index[population_idx], shell
                            ].to_numpy(dtype=np.float64)
                        )
                        estimated_population_sum = estimated_population.sum()
                        if (
                            np.isfinite(estimated_population).all()
                            and np.all(estimated_population >= 0.0)
                            and estimated_population_sum > 0.0
                        ):
                            populated_normalized_population[
                                populated_shell_idx
                            ] = estimated_population / estimated_population_sum
                            frozen_key = (int(atomic_number), shell)
                            if (
                                frozen_key
                                not in self._logged_frozen_ion_populations
                            ):
                                logger.warning(
                                    "Freezing ion populations for atomic "
                                    "number %s, shell %s because its rate "
                                    "matrix is singular; retaining the "
                                    "normalized pre-evaluator ion fractions.",
                                    atomic_number,
                                    shell,
                                )
                                self._logged_frozen_ion_populations.add(
                                    frozen_key
                                )
                            continue
                        if maximum_electron_density is not None:
                            maximum_shell_electron_density = (
                                maximum_electron_density[shell_idx]
                            )
                            maximum_charge_fraction = (
                                atomic_number
                                * elemental_density[shell_idx]
                                / maximum_shell_electron_density
                                if maximum_shell_electron_density > 0.0
                                else np.inf
                            )
                            if (
                                np.isfinite(maximum_charge_fraction)
                                and maximum_charge_fraction <= CHARGE_TOLERANCE
                            ):
                                skipped_populated_shells[
                                    populated_shell_idx
                                ] = True
                                populated_normalized_population[
                                    populated_shell_idx, 0
                                ] = 1.0
                                skipped_key = (int(atomic_number), shell)
                                if (
                                    skipped_key
                                    not in self._logged_negligible_ion_populations
                                ):
                                    logger.warning(
                                        "Skipping ion populations for atomic "
                                        "number %s, shell %s because its rate "
                                        "matrix is singular, it has no usable "
                                        "prior population, and its maximum "
                                        "possible charge fraction %.3e is "
                                        "below the charge tolerance %.3e.",
                                        atomic_number,
                                        shell,
                                        maximum_charge_fraction,
                                        CHARGE_TOLERANCE,
                                    )
                                    self._logged_negligible_ion_populations.add(
                                        skipped_key
                                    )
                                continue
                        singular_values = np.linalg.svd(
                            matrix, compute_uv=False
                        )
                        rank = np.linalg.matrix_rank(matrix)
                        condition_number = np.linalg.cond(matrix)
                        raise PlasmaIonizationError(
                            "Singular ion population matrix for atomic number "
                            f"{atomic_number}, shell {shell} "
                            f"(active_shell_idx={shell_idx}, "
                            "electron_density="
                            f"{electron_density[shell_idx]:.3e} cm^-3, "
                            f"rank={rank}/{matrix.shape[0]}, "
                            "smallest_singular_value="
                            f"{singular_values[-1]:.3e}, "
                            f"condition_number={condition_number:.3e})."
                        ) from exc
            nonfinite_populations = ~np.isfinite(
                populated_normalized_population
            ).all(axis=1)
            minimum_population = populated_normalized_population.min(axis=1)
            invalid_populations = nonfinite_populations | (
                minimum_population < -1e-12
            )
            if np.any(invalid_populations):
                populated_shell_idx = np.flatnonzero(invalid_populations)[0]
                shell_idx = populated_shell_indices[populated_shell_idx]
                raise PlasmaIonizationError(
                    "Nonfinite or negative ion population for atomic number "
                    f"{atomic_number}, shell {rate_matrices.columns[shell_idx]}."
                    f"{minimum_population[populated_shell_idx]}."
                )

            populated_normalized_population[
                populated_normalized_population < 0.0
            ] = 0.0
            populated_normalized_population /= (
                populated_normalized_population.sum(axis=1)[:, None]
            )
            populated_normalized_population[skipped_populated_shells] = 0.0
            normalized_population = np.zeros(
                (len(matrices), matrices.shape[1]), dtype=np.float64
            )
            normalized_population[populated_shell_indices] = (
                populated_normalized_population
            )

            ion_population[population_idx] = (
                normalized_population.T * elemental_density
            )

        self.rates_matrices = rate_matrices
        return ion_population

    def solve(
        self,
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField
        | None,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        elemental_number_density: pd.DataFrame,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame | float,
        boltzmann_factor: pd.DataFrame,
        tolerance: float = 1e-10,  # Started at 1e-14 - testing
        level_to_continuum_saha_factor: pd.DataFrame | None = None,
        lte_ionization_factor: pd.DataFrame | None = None,
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Solve ion populations while enforcing charge conservation."""
        if level_to_continuum_saha_factor is None:
            raise ValueError(
                "level_to_continuum_saha_factor is required when charge "
                "conservation is enabled."
            )
        bound_level_index = lte_level_population.index.get_level_values(
            "ion_number"
        ) < lte_level_population.index.get_level_values("atomic_number")
        return self._solve_charge_conserving(
            radiation_field,
            thermal_electron_energy_distribution,
            elemental_number_density,
            lte_level_population.loc[bound_level_index],
            estimated_level_population.loc[bound_level_index],
            lte_ion_population,
            estimated_ion_population,
            partition_function,
            boltzmann_factor.loc[bound_level_index],
            level_to_continuum_saha_factor.loc[
                lte_level_population.index[bound_level_index]
            ],
            tolerance,
            lte_ionization_factor=lte_ionization_factor,
        )

    def solve_charge_balance(
        self,
        electron_density: npt.NDArray[np.float64],
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField
        | None,
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
        lte_ionization_factor: pd.DataFrame | None = None,
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
        lte_ionization_factor : pandas.DataFrame, optional
            LTE ionization factors used to normalize charge rates.

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
            lte_ionization_factor=lte_ionization_factor,
            maximum_electron_density=maximum_electron_density,
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
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField
        | None,
        thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
        lte_level_population: pd.DataFrame,
        estimated_level_population: pd.DataFrame,
        lte_ion_population: pd.DataFrame,
        estimated_ion_population: pd.DataFrame,
        partition_function: pd.DataFrame,
        boltzmann_factor: pd.DataFrame,
        level_to_continuum_saha_factor: pd.DataFrame,
        elemental_number_density: pd.DataFrame,
        lte_ionization_factor: pd.DataFrame | None = None,
    ) -> float:
        """Solve the charge balance for one shell.

        Parameters
        ----------
        shell_idx : int
            Index of the shell whose electron density is being solved.
        maximum_electron_density : float
            Maximum possible electron number density in the shell.
        radiation_field : DilutePlanckianRadiationField | PlanckianRadiationField, optional
            Radiation field used to calculate ionization rates. Estimated rates
            do not require one.
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
        lte_ionization_factor : pandas.DataFrame, optional
            LTE ionization factors used to normalize charge rates.

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
        shell_columns = [elemental_number_density.columns[shell_idx]]
        if radiation_field is None:
            shell_radiation_field = None
        elif isinstance(radiation_field, DilutePlanckianRadiationField):
            shell_radiation_field = DilutePlanckianRadiationField(
                radiation_field.temperature[[shell_idx]],
                radiation_field.dilution_factor[[shell_idx]],
            )
        else:
            shell_radiation_field = PlanckianRadiationField(
                radiation_field.temperature[[shell_idx]]
            )
        shell_electron_distribution = ThermalElectronEnergyDistribution(
            thermal_electron_energy_distribution.energy,
            thermal_electron_energy_distribution.temperature[[shell_idx]],
            thermal_electron_energy_distribution.number_density[[shell_idx]],
        )
        shell_partition_function = (
            partition_function[shell_columns]
            if isinstance(partition_function, pd.DataFrame)
            else partition_function
        )
        shell_lte_level_population = lte_level_population[shell_columns]
        shell_estimated_level_population = estimated_level_population[
            shell_columns
        ]
        shell_lte_ion_population = lte_ion_population[shell_columns]
        shell_estimated_ion_population = estimated_ion_population[shell_columns]
        shell_boltzmann_factor = boltzmann_factor[shell_columns]
        shell_level_to_continuum_saha_factor = level_to_continuum_saha_factor[
            shell_columns
        ]
        shell_elemental_number_density = elemental_number_density[shell_columns]
        shell_maximum_electron_density = np.array([maximum_electron_density])
        shell_lte_ionization_factor = (
            None
            if lte_ionization_factor is None
            else lte_ionization_factor[shell_columns]
        )

        def charge_residual(
            electron_density_fraction: float,
        ) -> float:
            """Calculate the normalized charge residual for one trial density."""
            electron_density = np.array(
                [electron_density_fraction * maximum_electron_density]
            )
            return self.solve_charge_balance(
                electron_density,
                shell_radiation_field,
                shell_electron_distribution,
                shell_lte_level_population,
                shell_estimated_level_population,
                shell_lte_ion_population,
                shell_estimated_ion_population,
                shell_partition_function,
                shell_boltzmann_factor,
                shell_level_to_continuum_saha_factor,
                shell_elemental_number_density,
                shell_maximum_electron_density,
                lte_ionization_factor=shell_lte_ionization_factor,
            )[1][0]

        try:
            electron_density_fraction, root_result = brentq(
                charge_residual,
                MINIMUM_ELECTRON_DENSITY_FRACTION,
                1.0,
                xtol=CHARGE_TOLERANCE,
                rtol=MINIMUM_BRENT_RELATIVE_TOLERANCE,
                maxiter=self.max_solver_iterations,
                full_output=True,
            )
        except ValueError as exc:
            lower_residual = charge_residual(MINIMUM_ELECTRON_DENSITY_FRACTION)
            upper_residual = charge_residual(1.0)
            raise PlasmaIonizationError(
                f"Charge residual does not bracket shell {shell_idx}: "
                "Q_hat(near-neutral)="
                f"{lower_residual}, Q_hat(1)={upper_residual}."
            ) from exc

        shell = elemental_number_density.columns[shell_idx]
        terminal_electron_density = np.array(
            [electron_density_fraction * maximum_electron_density]
        )
        _, terminal_charge_residual = self.solve_charge_balance(
            terminal_electron_density,
            shell_radiation_field,
            shell_electron_distribution,
            shell_lte_level_population,
            shell_estimated_level_population,
            shell_lte_ion_population,
            shell_estimated_ion_population,
            shell_partition_function,
            shell_boltzmann_factor,
            shell_level_to_continuum_saha_factor,
            shell_elemental_number_density,
            shell_maximum_electron_density,
            lte_ionization_factor=shell_lte_ionization_factor,
        )
        self._charge_root_diagnostics[shell] = (
            terminal_charge_residual[0],
            int(root_result.iterations),
            int(root_result.function_calls),
            bool(root_result.converged),
            str(root_result.flag),
        )
        return terminal_electron_density[0]

    def _solve_charge_conserving(
        self,
        radiation_field: DilutePlanckianRadiationField
        | PlanckianRadiationField
        | None,
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
        lte_ionization_factor: pd.DataFrame | None = None,
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
        lte_ionization_factor : pandas.DataFrame, optional
            LTE ionization factors used to normalize charge rates.

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
            logger.debug("Ion solver iteration %d", iteration + 1)
            self._charge_root_diagnostics = {}
            prepare_rate_matrix = getattr(
                self.rate_matrix_solver, "prepare", None
            )
            if prepare_rate_matrix is not None:
                rate_matrix_mode = "prepared"
                logger.debug(
                    "Ion solver iteration %d starting ion rate-matrix "
                    "preparation.",
                    iteration + 1,
                )
                preparation_start_time = perf_counter()
                prepare_rate_matrix(
                    thermal_electron_energy_distribution,
                    estimated_level_population,
                    estimated_ion_population,
                    partition_function,
                    boltzmann_factor,
                    level_to_continuum_saha_factor,
                    lte_ionization_factor,
                )
                preparation_elapsed = perf_counter() - preparation_start_time
                logger.log(
                    logging.INFO
                    if preparation_elapsed >= 10.0
                    else logging.DEBUG,
                    "Ion solver iteration %d finished ion rate-matrix "
                    "preparation in %.3f s.",
                    iteration + 1,
                    preparation_elapsed,
                )
            else:
                rate_matrix_mode = "direct"
                logger.debug(
                    "Ion solver iteration %d uses direct ion rate-matrix "
                    "construction (%s).",
                    iteration + 1,
                    type(self.rate_matrix_solver).__name__,
                )
            unconverged_shell_indices = np.flatnonzero(~converged_shells)
            logger.debug(
                "Ion solver iteration %d starting %d shell charge solves "
                "(rate-matrix mode=%s).",
                iteration + 1,
                len(unconverged_shell_indices),
                rate_matrix_mode,
            )
            shell_charge_start_time = perf_counter()
            for shell_idx in unconverged_shell_indices:
                electron_density[shell_idx] = self.solve_shell_charge(
                    shell_idx,
                    maximum_electron_density_array[shell_idx],
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
                    lte_ionization_factor=lte_ionization_factor,
                )
            charge_solve_elapsed = perf_counter() - shell_charge_start_time
            logger.log(
                logging.INFO if charge_solve_elapsed >= 10.0 else logging.DEBUG,
                "Ion solver iteration %d finished %d shell charge solves in "
                "%.3f s (rate-matrix mode=%s).",
                iteration + 1,
                len(unconverged_shell_indices),
                charge_solve_elapsed,
                rate_matrix_mode,
            )
            logger.debug(
                "Ion solver iteration %d starting full charge-balance "
                "reconstruction.",
                iteration + 1,
            )
            charge_balance_start_time = perf_counter()
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
                    lte_ionization_factor=lte_ionization_factor,
                )
            )
            logger.debug(
                "Ion solver iteration %d finished full charge-balance "
                "reconstruction in %.3f s.",
                iteration + 1,
                perf_counter() - charge_balance_start_time,
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
            active_shell_indices = np.flatnonzero(~converged_shells)
            if len(estimated_population_indices) == 0:
                population_converged = np.ones_like(
                    converged_shells, dtype=bool
                )
                maximum_population_residual = 0.0
                maximum_population_ion_stage = None
                maximum_population_shell = None
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
                active_population_residual = np.abs(
                    population_delta[:, active_shell_indices]
                )
                ranking_residual = np.nan_to_num(
                    active_population_residual,
                    nan=np.inf,
                    posinf=np.inf,
                    neginf=np.inf,
                )
                (
                    population_position,
                    active_shell_position,
                ) = np.unravel_index(
                    np.argmax(ranking_residual), ranking_residual.shape
                )
                shell_idx = active_shell_indices[active_shell_position]
                estimated_population_idx = estimated_population_indices[
                    population_position
                ]
                maximum_population_residual = active_population_residual[
                    population_position, active_shell_position
                ]
                maximum_population_ion_stage = tuple(
                    int(index_value)
                    for index_value in estimated_ion_population.index[
                        estimated_population_idx
                    ]
                )
                maximum_population_shell = elemental_number_density.columns[
                    shell_idx
                ]

            active_charge_residual = np.abs(
                charge_residual[active_shell_indices]
            )
            charge_position = int(
                np.argmax(
                    np.nan_to_num(
                        active_charge_residual,
                        nan=np.inf,
                        posinf=np.inf,
                        neginf=np.inf,
                    )
                )
            )
            maximum_charge_residual = active_charge_residual[charge_position]
            maximum_charge_shell = elemental_number_density.columns[
                active_shell_indices[charge_position]
            ]
            charge_converged = np.abs(charge_residual) < CHARGE_TOLERANCE
            population_failure_mask = ~population_converged[
                active_shell_indices
            ]
            charge_failure_mask = ~charge_converged[active_shell_indices]
            population_failure_count = np.count_nonzero(population_failure_mask)
            charge_failure_count = np.count_nonzero(charge_failure_mask)
            joint_failure_count = np.count_nonzero(
                population_failure_mask & charge_failure_mask
            )
            population_blocked = population_failure_count != 0
            charge_blocked = charge_failure_count != 0
            if population_blocked and charge_blocked:
                blocker = "population and charge"
            elif population_blocked:
                blocker = "population"
            elif charge_blocked:
                blocker = "charge"
            else:
                blocker = "none"

            converged_shells |= charge_converged & population_converged
            logger.debug(
                "Ion solver iteration %d: %d/%d shells unconverged; "
                "population-failing shells=%d, charge-failing shells=%d, "
                "failures of both=%d; maximum ion-stage relative residual "
                "%.3e at ion stage %s in shell %s (tolerance=%.3e); maximum "
                "normalized charge residual %.3e in shell %s "
                "(tolerance=%.3e); blocker=%s",
                iteration + 1,
                np.count_nonzero(~converged_shells),
                len(converged_shells),
                population_failure_count,
                charge_failure_count,
                joint_failure_count,
                maximum_population_residual,
                maximum_population_ion_stage,
                maximum_population_shell,
                tolerance,
                maximum_charge_residual,
                maximum_charge_shell,
                CHARGE_TOLERANCE,
                blocker,
            )
            charge_root_diagnostic = self._charge_root_diagnostics.get(
                maximum_charge_shell
            )
            if charge_root_diagnostic is None:
                logger.warning(
                    "Ion solver iteration %d has no charge-root diagnostic "
                    "for shell %s.",
                    iteration + 1,
                    maximum_charge_shell,
                )
            else:
                (
                    local_charge_residual,
                    root_iterations,
                    root_function_calls,
                    root_converged,
                    root_flag,
                ) = charge_root_diagnostic
                final_charge_residual = charge_residual[
                    active_shell_indices[charge_position]
                ]
                logger.debug(
                    "Ion solver iteration %d charge-root consistency in shell "
                    "%s: local terminal residual=%.3e, final residual=%.3e, "
                    "absolute difference=%.3e; Brent converged=%s, "
                    "iterations=%d, function_calls=%d, flag=%s",
                    iteration + 1,
                    maximum_charge_shell,
                    local_charge_residual,
                    final_charge_residual,
                    abs(local_charge_residual - final_charge_residual),
                    root_converged,
                    root_iterations,
                    root_function_calls,
                    root_flag,
                )
            if np.all(converged_shells):
                logger.debug(
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
                estimated_ion_population = estimated_ion_population.copy()
                estimated_ion_population.iloc[estimated_population_indices] = (
                    ion_population_solution[solution_population_indices]
                )
                level_ion_index = estimated_level_population.index.droplevel(
                    "level_number"
                )
                level_solution_indices = (
                    self.rate_matrix_solver.ion_population_index.get_indexer(
                        level_ion_index
                    )
                )
                solved_levels = level_solution_indices >= 0
                level_totals = estimated_level_population.groupby(
                    level=["atomic_number", "ion_number"]
                ).transform("sum")
                level_fractions = estimated_level_population.divide(
                    level_totals.where(level_totals != 0.0, 1.0)
                )
                estimated_level_population = estimated_level_population.copy()
                estimated_level_population.iloc[solved_levels] = (
                    level_fractions.to_numpy()[solved_levels]
                    * ion_population_solution[
                        level_solution_indices[solved_levels]
                    ]
                )

        raise PlasmaIonizationError(
            "Ion population solver did not converge after "
            f"{self.max_solver_iterations} iterations."
        )


class FixedElectronDensityIonPopulationSolver(IonPopulationSolver):
    """Solve ion populations for a fixed electron density."""

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
        tolerance: float = 1e-14,
        lte_ionization_factor: pd.DataFrame | None = None,
    ) -> tuple[pd.DataFrame, pd.Series]:
        """Solve ion populations without imposing charge conservation."""
        bound_level_index = lte_level_population.index.get_level_values(
            "ion_number"
        ) < lte_level_population.index.get_level_values("atomic_number")
        return self._solve_fixed_electron_density(
            radiation_field,
            thermal_electron_energy_distribution,
            elemental_number_density,
            lte_level_population.loc[bound_level_index],
            estimated_level_population.loc[bound_level_index],
            lte_ion_population,
            estimated_ion_population,
            partition_function,
            boltzmann_factor.loc[bound_level_index],
            tolerance,
            lte_ionization_factor=lte_ionization_factor,
        )

    def _solve_fixed_electron_density(
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
        lte_ionization_factor: pd.DataFrame | None = None,
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
        lte_ionization_factor : pandas.DataFrame, optional
            LTE ionization factors used to normalize ionization rates.

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
                lte_ionization_factor=lte_ionization_factor,
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
                np.vstack(solved_matrices.to_numpy()[0]).T,
                index=self.rate_matrix_solver.ion_population_index,
                columns=self.rates_matrices.columns,
            )

            if (ion_population_solution < 0).any().any():
                ion_population_solution[ion_population_solution < 0] = 0.0

            electron_population_solution = (
                ion_population_solution
                * ion_population_solution.index.get_level_values("ion_number")
                .to_numpy()[np.newaxis]
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
                electron_population_solution.to_numpy() * u.cm**-3
            )
        else:
            logger.warning(
                "Ion population solver did not converge after %d iterations.",
                iteration,
            )

        return ion_population_solution, electron_population_solution
