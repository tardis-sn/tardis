import logging

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.optimize import OptimizeResult, least_squares
from scipy.sparse import block_diag

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.evaluator import (
    PlasmaEquilibriumEvaluation,
    PlasmaEquilibriumEvaluator,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField

logger = logging.getLogger(__name__)

THERMAL_BALANCE_RESIDUAL_LIMITS = {
    "trial_level_residual": 1e-10,
    "level_residual": 1e-10,
    "charge_residual": 1e-10,
    # The charge projection independently solves for the electron density, so
    # the projected value can differ slightly from the outer trial value.
    "electron_residual": 2e-8,
    # Heating residuals subtract nearly cancelling rates, leaving a small
    # floating-point floor.
    "total_heating": 5e-13,
    "fractional_heating": 2e-7,
}


class ThermalBalanceSolver:
    """Class to solve the thermal balance equation using all available

    heating and cooling processes. See section 6.4 in Lucy 03, especially
    equations 55-57.
    """

    def __init__(
        self,
        bound_free_solver: BoundFreeThermalRates,
        free_free_solver: FreeFreeThermalRates,
        collisional_ionization_solver: CollisionalIonizationThermalRates,
        collisional_bound_solver: CollisionalBoundThermalRates,
    ) -> None:
        """Initialize the thermal balance solver with individual solver components.

        Parameters
        ----------
        bound_free_solver
            Solver for bound-free heating and free-bound cooling processes.
        free_free_solver
            Solver for free-free heating and cooling processes.
        collisional_ionization_solver
            Solver for collisional ionization heating and cooling processes.
        collisional_bound_solver
            Solver for collisional bound-bound heating and cooling processes.
        """
        self.bound_free_solver = bound_free_solver
        self.free_free_solver = free_free_solver
        self.collisional_ionization_solver = collisional_ionization_solver
        self.collisional_bound_solver = collisional_bound_solver

    def solve(
        self,
        thermal_electron_distribution: ThermalElectronEnergyDistribution,
        level_population: pd.DataFrame,
        ion_population: pd.DataFrame,
        collisional_ionization_rate_coefficient: pd.DataFrame,
        collisional_deexcitation_rate_coefficient: pd.DataFrame,
        collisional_excitation_rate_coefficient: pd.DataFrame,
        free_free_heating_estimator: pd.DataFrame,
        level_population_ratio: pd.DataFrame,
        radiation_field: DilutePlanckianRadiationField | None = None,
        bound_free_heating_estimator: pd.DataFrame | None = None,
        stimulated_recombination_estimator: pd.DataFrame | None = None,
    ) -> tuple[pd.Series, pd.Series]:
        """Compute the current heating rate and the fractional heating

        rate using all available processes. See section 6.4 in Lucy 03.

        Parameters
        ----------
        thermal_electron_distribution : ThermalElectronEnergyDistribution
            Electron energy, temperature, and density.
        level_population : pd.DataFrame
            Level number density.
        ion_population : pd.DataFrame
            Ion number density.
        collisional_ionization_rate_coefficient : pd.DataFrame
            Collisional ionization rate coefficient.
        collisional_deexcitation_rate_coefficient : pd.DataFrame
            Collisional deexcitation rate coefficient.
        collisional_excitation_rate_coefficient : pd.DataFrame
            Collisional excitation rate coefficient.
        free_free_heating_estimator : pd.DataFrame
            Montecarlo estimator for free-free heating.
        level_population_ratio : pd.DataFrame
            Level population to ion population ratio. Lucy 03, equation 14.
        radiation_field : RadiationField, optional
            Radiation field for mean intensity calculation.
        bound_free_heating_estimator : pd.DataFrame, optional
            Bound-free heating estimator.
        stimulated_recombination_estimator : pd.DataFrame, optional
            Stimulated recombination estimator.

        Returns
        -------
        tuple[pd.Series, pd.Series]
            Total heating rate and fractional heating rate for each cell.
        """
        electron_density = thermal_electron_distribution.number_density

        bound_free_heating, free_bound_cooling = self.bound_free_solver.solve(
            level_population,
            ion_population,
            thermal_electron_distribution,
            level_population_ratio,
            radiation_field,
            bound_free_heating_estimator,
            stimulated_recombination_estimator,
        )

        free_free_heating, free_free_cooling = self.free_free_solver.solve(
            free_free_heating_estimator,
            thermal_electron_distribution,
            ion_population,
        )

        collisional_ionization_heating, collisional_ionization_cooling = (
            self.collisional_ionization_solver.solve(
                electron_density,
                ion_population,
                level_population,
                collisional_ionization_rate_coefficient,
                level_population_ratio,
            )
        )

        collisional_bound_heating, collisional_bound_cooling = (
            self.collisional_bound_solver.solve(
                electron_density,
                collisional_deexcitation_rate_coefficient,
                collisional_excitation_rate_coefficient,
                level_population,
            )
        )

        total_heating = (
            bound_free_heating.to_numpy()
            + free_free_heating.to_numpy()
            + collisional_ionization_heating.to_numpy()
            + collisional_bound_heating.to_numpy()
        )
        total_cooling = (
            free_bound_cooling.to_numpy()
            + free_free_cooling.to_numpy()
            + collisional_ionization_cooling.to_numpy()
            + collisional_bound_cooling.to_numpy()
        )

        total_heating_rate = total_heating - total_cooling
        fractional_heating_rate = (
            total_heating - total_cooling
        ) / total_cooling

        return (
            pd.Series(total_heating_rate, index=level_population.columns),
            pd.Series(fractional_heating_rate, index=level_population.columns),
        )


class ElectronThermalBalanceSolver:
    """Solve electron density and temperature balance across plasma shells."""

    least_squares_settings = {
        "xtol": 1e-14,
        "ftol": 1e-12,
        "gtol": 1e-14,
        "x_scale": "jac",
        "max_nfev": 100,
    }

    def __init__(
        self,
        evaluator: PlasmaEquilibriumEvaluator,
        radiation_temperature: npt.NDArray[np.float64],
        level_initial_guess: pd.DataFrame,
        maximum_electron_density: npt.NDArray[np.float64],
    ) -> None:
        """Initialize a solve from one frozen plasma-equilibrium snapshot.

        Parameters
        ----------
        evaluator : PlasmaEquilibriumEvaluator
            Fixed equilibrium evaluator for the outer thermal solve.
        radiation_temperature : numpy.ndarray
            Radiation temperature in each shell [K].
        level_initial_guess : pandas.DataFrame
            Frozen normalized level populations used during optimization.
        maximum_electron_density : numpy.ndarray
            Maximum electron number density in each shell [cm^-3].
        """
        self.evaluator = evaluator
        self.radiation_temperature = radiation_temperature
        self.level_initial_guess = level_initial_guess
        self.maximum_electron_density = maximum_electron_density

    def calculate_residual(
        self, candidate: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """
        Calculate charge and heating residuals for all shells.

        Parameters
        ----------
        candidate : numpy.ndarray
            Interleaved electron-density fractions and electron-to-radiation
            temperature ratios for each shell.

        Returns
        -------
        numpy.ndarray
            Interleaved normalized electron-density and fractional-heating
            residuals for each shell.
        """
        electron_density_fraction = candidate[::2]
        temperature_ratio = candidate[1::2]

        logger.info("Link: %s", temperature_ratio)
        evaluation = self.evaluator.evaluate(
            self.maximum_electron_density * electron_density_fraction,
            self.radiation_temperature * temperature_ratio,
            self.level_initial_guess,
        )

        electron_residual = evaluation.electron_residual.to_numpy()
        fractional_heating = evaluation.fractional_heating.to_numpy()
        residual = np.zeros_like(candidate)
        if not np.isfinite(fractional_heating).all():
            logger.warning("Heating not finite\n")
        if not np.isfinite(electron_residual).all():
            logger.warning("Fractional electron change not finite\n")

        residual[::2] = electron_residual
        residual[1::2] = fractional_heating
        logger.info(
            "Normalized electron fraction change: %s", electron_residual
        )
        logger.info("Heating: %s", fractional_heating)
        return residual

    def solve_multiple_shells(
        self,
        initial_guess: npt.NDArray[np.float64],
        lower_bound: npt.NDArray[np.float64],
        upper_bound: npt.NDArray[np.float64],
    ) -> OptimizeResult:
        """
        Solve charge and heating balance for all shells jointly.

        Parameters
        ----------
        initial_guess : numpy.ndarray
            Initial interleaved electron-density fractions and temperature
            ratios for each shell.
        lower_bound : numpy.ndarray
            Lower bounds for each value in `initial_guess`.
        upper_bound : numpy.ndarray
            Upper bounds for each value in `initial_guess`.

        Returns
        -------
        scipy.optimize.OptimizeResult
            Result of the joint nonlinear least-squares solve.
        """
        no_shells = len(initial_guess) // 2
        return least_squares(
            self.calculate_residual,
            initial_guess,
            bounds=(lower_bound, upper_bound),
            jac_sparsity=block_diag([np.ones((2, 2))] * no_shells),
            verbose=1,
            method="trf",
            **self.least_squares_settings,
        )

    def solve_single_shell(
        self,
        shell_idx: int,
        accepted_candidate: npt.NDArray[np.float64],
        initial_guess: npt.NDArray[np.float64],
        lower_bound: npt.NDArray[np.float64],
        upper_bound: npt.NDArray[np.float64],
    ) -> OptimizeResult:
        """
        Retry one failed shell, including its alternate initial state.

        Parameters
        ----------
        shell_idx : int
            Zero-based position of the shell to solve.
        accepted_candidate : numpy.ndarray
            Current interleaved candidate for all shells. Values outside the
            selected shell remain fixed during the solve.
        initial_guess : numpy.ndarray
            Original bounded candidate used as the alternate starting state.
        lower_bound : numpy.ndarray
            Lower bounds for the interleaved candidate values.
        upper_bound : numpy.ndarray
            Upper bounds for the interleaved candidate values.

        Returns
        -------
        scipy.optimize.OptimizeResult
            Preferred result from the primary or alternate single-shell solve.
        """
        shell_slice = slice(2 * shell_idx, 2 * shell_idx + 2)
        local_lower_bound = lower_bound[shell_slice].copy()
        local_upper_bound = upper_bound[shell_slice]
        local_residual_limits = np.array(
            [
                THERMAL_BALANCE_RESIDUAL_LIMITS["electron_residual"],
                THERMAL_BALANCE_RESIDUAL_LIMITS["fractional_heating"],
            ]
        )

        def shell_residual(
            local_candidate: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            trial_candidate = accepted_candidate.copy()
            trial_candidate[shell_slice] = local_candidate
            return self.calculate_residual(trial_candidate)[shell_slice]

        local_start = np.clip(
            accepted_candidate[shell_slice],
            local_lower_bound,
            local_upper_bound,
        )
        shell_result = least_squares(
            shell_residual,
            local_start,
            bounds=(local_lower_bound, local_upper_bound),
            **self.least_squares_settings,
        )
        if np.any(np.abs(shell_result.fun) > local_residual_limits):
            alternate_start = np.clip(
                initial_guess[shell_slice],
                local_lower_bound,
                local_upper_bound,
            )
            if not np.array_equal(alternate_start, local_start):
                alternate_result = least_squares(
                    shell_residual,
                    alternate_start,
                    bounds=(local_lower_bound, local_upper_bound),
                    **self.least_squares_settings,
                )
                shell_closes = np.all(
                    np.abs(shell_result.fun) <= local_residual_limits
                )
                alternate_closes = np.all(
                    np.abs(alternate_result.fun) <= local_residual_limits
                )
                if (alternate_closes and not shell_closes) or (
                    alternate_closes == shell_closes
                    and np.linalg.norm(
                        alternate_result.fun / local_residual_limits
                    )
                    < np.linalg.norm(shell_result.fun / local_residual_limits)
                ):
                    shell_result = alternate_result

        logger.info(
            "Shell-local thermal balance shell=%s candidate=%s residual=%s status=%s",
            shell_idx,
            shell_result.x,
            shell_result.fun,
            shell_result.status,
        )
        return shell_result

    def evaluate_candidate(
        self, candidate: npt.NDArray[np.float64]
    ) -> PlasmaEquilibriumEvaluation:
        """
        Rebuild the terminal equilibrium state for a candidate.

        Parameters
        ----------
        candidate : numpy.ndarray
            Interleaved electron-density fractions and temperature ratios for
            each shell.

        Returns
        -------
        tardis.plasma.equilibrium.evaluator.PlasmaEquilibriumEvaluation
            Equilibrium state rebuilt using the candidate's converged level
            populations as the terminal seed.
        """
        electron_density = self.maximum_electron_density * candidate[::2]
        electron_temperature = np.asarray(
            self.radiation_temperature * candidate[1::2], dtype=np.float64
        )
        seed_evaluation = self.evaluator.evaluate(
            electron_density,
            electron_temperature,
            self.level_initial_guess,
        )
        return self.evaluator.evaluate(
            electron_density,
            electron_temperature,
            seed_evaluation.normalized_population,
        )

    @staticmethod
    def validate_evaluation(evaluation: PlasmaEquilibriumEvaluation) -> None:
        """
        Validate population, charge, level, and thermal closure.

        Parameters
        ----------
        evaluation : tardis.plasma.equilibrium.evaluator.PlasmaEquilibriumEvaluation
            Terminal equilibrium state to validate.

        Raises
        ------
        RuntimeError
            If populations or derived states are invalid, or if any closure
            residual exceeds its declared tolerance.
        """
        failures = []
        population_fields = {
            "normalized_population": evaluation.normalized_population,
            "absolute_level_population": evaluation.absolute_level_population,
            "ion_population": evaluation.ion_population,
            "electron_density": evaluation.charge_solved_electron_density,
        }
        for field_name, population in population_fields.items():
            if population is None:
                failures.append(f"{field_name} is missing")
                continue
            population_values = np.asarray(population, dtype=np.float64)
            if not np.isfinite(population_values).all():
                failures.append(f"{field_name} is nonfinite")
            if (population_values < 0.0).any():
                failures.append(f"{field_name} is negative")

        normalized_totals = evaluation.normalized_population.sum(axis=0)
        if not np.allclose(normalized_totals, 1.0, rtol=1e-12, atol=0.0):
            failures.append("normalized_population does not sum to one")

        for field_name, state in {
            "tau_sobolev": evaluation.tau_sobolev,
            "beta_sobolev": evaluation.beta_sobolev,
        }.items():
            if not np.isfinite(np.asarray(state, dtype=np.float64)).all():
                failures.append(f"{field_name} is nonfinite")

        for field_name, tolerance in THERMAL_BALANCE_RESIDUAL_LIMITS.items():
            residual = getattr(evaluation, field_name)
            if residual is None:
                failures.append(f"{field_name} is missing")
                continue
            residual_values = np.asarray(residual, dtype=np.float64)
            if not np.isfinite(residual_values).all():
                failures.append(f"{field_name} is nonfinite")
                continue
            maximum_residual = float(np.max(np.abs(residual_values)))
            if maximum_residual > tolerance:
                failures.append(
                    f"{field_name}={maximum_residual:.3e} exceeds "
                    f"{tolerance:.3e}"
                )

        if failures:
            raise RuntimeError(
                "Thermal-balance solution failed final closure: "
                + "; ".join(failures)
            )

    def solve(
        self,
        initial_guess: npt.NDArray[np.float64],
        lower_bound: npt.NDArray[np.float64],
        upper_bound: npt.NDArray[np.float64],
    ) -> tuple[npt.NDArray[np.float64], PlasmaEquilibriumEvaluation]:
        """
        Solve jointly, then retry shells that fail final closure.

        Parameters
        ----------
        initial_guess : numpy.ndarray
            Initial interleaved electron-density fractions and temperature
            ratios for each shell.
        lower_bound : numpy.ndarray
            Lower bounds for each value in `initial_guess`.
        upper_bound : numpy.ndarray
            Upper bounds for each value in `initial_guess`.

        Returns
        -------
        accepted_candidate : numpy.ndarray
            Interleaved electron-density fractions and temperature ratios from
            the accepted joint and single-shell solves.
        evaluation : tardis.plasma.equilibrium.evaluator.PlasmaEquilibriumEvaluation
            Validated terminal equilibrium state for `accepted_candidate`.

        Raises
        ------
        RuntimeError
            If the joint failure cannot be attributed to individual charge or
            heating residuals, or if shell-local retries do not restore final
            closure.
        """
        multiple_shell_result = self.solve_multiple_shells(
            initial_guess, lower_bound, upper_bound
        )
        accepted_candidate = multiple_shell_result.x
        evaluation = self.evaluate_candidate(accepted_candidate)
        try:
            self.validate_evaluation(evaluation)
        except RuntimeError as joint_failure:
            logger.warning(
                "Joint thermal balance failed final closure; retrying failed "
                "electron/heating shells: %s",
                joint_failure,
            )
            local_residual_limits = np.array(
                [
                    THERMAL_BALANCE_RESIDUAL_LIMITS["electron_residual"],
                    THERMAL_BALANCE_RESIDUAL_LIMITS["fractional_heating"],
                ]
            )
            electron_residual = evaluation.electron_residual.to_numpy()
            fractional_heating = evaluation.fractional_heating.to_numpy()
            failed_shells = np.flatnonzero(
                ~np.isfinite(electron_residual)
                | ~np.isfinite(fractional_heating)
                | (np.abs(electron_residual) > local_residual_limits[0])
                | (np.abs(fractional_heating) > local_residual_limits[1])
            )
            if not len(failed_shells):
                raise
            accepted_candidate = accepted_candidate.copy()
            for shell_idx in failed_shells:
                shell_result = self.solve_single_shell(
                    shell_idx,
                    accepted_candidate,
                    initial_guess,
                    lower_bound,
                    upper_bound,
                )
                shell_slice = slice(2 * shell_idx, 2 * shell_idx + 2)
                accepted_candidate[shell_slice] = shell_result.x

            evaluation = self.evaluate_candidate(accepted_candidate)
            self.validate_evaluation(evaluation)

        return accepted_candidate, evaluation
