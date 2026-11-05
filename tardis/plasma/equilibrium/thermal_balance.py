import logging
from collections.abc import Callable, Iterable

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.optimize import OptimizeResult, least_squares
from scipy.sparse import block_diag

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.evaluator import (
    LEVEL_RESIDUAL_LIMIT,
    PlasmaEquilibriumEvaluation,
    PlasmaEquilibriumEvaluator,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.exceptions import PlasmaIonizationError
from tardis.plasma.radiation_field import DilutePlanckianRadiationField

logger = logging.getLogger(__name__)

# Inner closure required of every accepted thermal state. The electron and
# heating residuals close at their Monte Carlo noise instead.
THERMAL_BALANCE_RESIDUAL_LIMITS = {
    "relative_trial_level_residual": LEVEL_RESIDUAL_LIMIT,
    "charge_residual": 1e-10,
}

# Relative step of the forward differences that measure the Jacobian
# diagonal used to scale each residual row. Only the magnitude of the
# derivative matters, and the residual is smooth to ~1e-9 at this step.
ROW_SCALE_RELATIVE_STEP = 1e-4

# A noise-limited closure above this level means the plasma state is too
# poorly constrained by the Monte Carlo estimators.
NOISE_LIMITED_CLOSURE_WARNING = 1e-2


def calculate_batch_standard_error(
    batch_values: npt.ArrayLike,
    batch_packet_counts: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """Calculate the standard error of a pooled Monte Carlo estimate.

    Each batch value is the estimate from its own independent packets. The
    pooled estimate is the packet-weighted mean of the batches, and its
    variance is estimated from the scatter of the batches as

    ``sum_b n_b (x_b - x_pooled)**2 / (N (K - 1))``,

    which reduces to the variance of the batch means divided by ``K`` for
    equal batches.

    Parameters
    ----------
    batch_values : array-like
        Batch estimates with the batch along the first axis.
    batch_packet_counts : array-like
        Number of packets in each batch.

    Returns
    -------
    numpy.ndarray
        Standard error of the pooled estimate, with the shape of one batch.

    Raises
    ------
    ValueError
        If fewer than two batches are given, or the number of packet counts
        differs from the number of batches.
    """
    batch_values = np.asarray(batch_values, dtype=np.float64)
    batch_packet_counts = np.asarray(batch_packet_counts, dtype=np.float64)
    number_of_batches = len(batch_packet_counts)
    if number_of_batches < 2:
        raise ValueError("At least two batches are required.")
    if len(batch_values) != number_of_batches:
        raise ValueError("Each batch needs exactly one packet count.")
    total_packet_count = batch_packet_counts.sum()
    pooled_value = np.tensordot(
        batch_packet_counts / total_packet_count, batch_values, axes=1
    )
    weighted_scatter = np.tensordot(
        batch_packet_counts, (batch_values - pooled_value) ** 2, axes=1
    )
    return np.sqrt(
        weighted_scatter / (total_packet_count * (number_of_batches - 1))
    )


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
    """Solve electron density and temperature balance across plasma shells.

    The electron and heating residuals of each shell close at twice their
    Monte Carlo standard error, measured from independent packet batches.
    """

    # max_nfev counts the optimizer's own residual evaluations; its
    # finite-difference Jacobian columns are evaluated in addition.
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
        batch_evaluators: Callable[[], Iterable[PlasmaEquilibriumEvaluator]],
        batch_packet_counts: npt.NDArray[np.int64],
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
        batch_evaluators : callable
            Returns evaluators of independent Monte Carlo packet batches,
            whose scatter sets the closure limits.
        batch_packet_counts : numpy.ndarray
            Number of packets in each batch of ``batch_evaluators``.

        Raises
        ------
        ValueError
            If fewer than two packet batches are given.
        """
        if len(batch_packet_counts) < 2:
            raise ValueError(
                "Thermal balance needs at least two Monte Carlo packet "
                "batches to measure its closure limits."
            )
        self.evaluator = evaluator
        self.radiation_temperature = radiation_temperature
        self.level_initial_guess = level_initial_guess
        self.maximum_electron_density = maximum_electron_density
        self.batch_evaluators = batch_evaluators
        self.batch_packet_counts = batch_packet_counts
        # Set by ``solve``: the 2-sigma limits of the outer residuals.
        self.outer_residual_limits: npt.NDArray[np.float64] | None = None

    @staticmethod
    def calculate_electron_and_thermal_residuals(
        evaluation: PlasmaEquilibriumEvaluation,
    ) -> npt.NDArray[np.float64]:
        """Interleave the residuals that the thermal solve drives to zero.

        The outer residuals of each shell are the normalized electron-density
        change (n_e,solved - n_e,trial)/n_e,trial and the fractional heating
        (H - C)/C. Level and charge residuals are solved inside each
        evaluation.

        Parameters
        ----------
        evaluation : tardis.plasma.equilibrium.evaluator.PlasmaEquilibriumEvaluation
            Equilibrium state at a candidate.

        Returns
        -------
        numpy.ndarray
            Interleaved electron and heating residuals for each shell.
        """
        outer_residuals = np.empty(2 * len(evaluation.electron_residual))
        outer_residuals[::2] = evaluation.electron_residual.to_numpy()
        outer_residuals[1::2] = evaluation.fractional_heating.to_numpy()
        return outer_residuals

    def calculate_residual(
        self,
        candidate: npt.NDArray[np.float64],
        evaluator: PlasmaEquilibriumEvaluator | None = None,
    ) -> npt.NDArray[np.float64]:
        """
        Calculate the outer residuals for all shells at a candidate.

        The equilibrium is evaluated from the frozen level seed.

        Parameters
        ----------
        candidate : numpy.ndarray
            Interleaved electron-density fractions and electron-to-radiation
            temperature ratios for each shell.
        evaluator : PlasmaEquilibriumEvaluator, optional
            Evaluator to use instead of the solver's own.

        Returns
        -------
        numpy.ndarray
            Interleaved electron and heating residuals for each shell (see
            ``calculate_outer_residuals``).
        """
        evaluator = self.evaluator if evaluator is None else evaluator
        logger.debug("Link: %s", candidate[1::2])
        evaluation = evaluator.evaluate(
            self.maximum_electron_density * candidate[::2],
            self.radiation_temperature * candidate[1::2],
            self.level_initial_guess,
        )
        residual = self.calculate_electron_and_thermal_residuals(evaluation)
        if not np.isfinite(residual[1::2]).all():
            logger.warning("Heating not finite\n")
        if not np.isfinite(residual[::2]).all():
            logger.warning("Fractional electron change not finite\n")
        logger.debug("Normalized electron fraction change: %s", residual[::2])
        logger.debug("Heating: %s", residual[1::2])
        return residual

    def calculate_trial_residual(
        self, candidate: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Calculate residuals at an optimizer trial point.

        A trial point at which the ionization balance cannot be solved is
        rejected with nonfinite residuals. At a trust-region trial step this
        makes ``least_squares`` shorten its step. At the initial point the
        optimizer raises, and a nonfinite finite-difference column makes its
        Jacobian nonfinite.

        Parameters
        ----------
        candidate : numpy.ndarray
            Interleaved electron-density fractions and electron-to-radiation
            temperature ratios for each shell.

        Returns
        -------
        numpy.ndarray
            Residuals from ``calculate_residual``, or NaN if the ionization
            balance failed at the candidate.
        """
        try:
            return self.calculate_residual(candidate)
        except PlasmaIonizationError as error:
            logger.warning(
                "Rejecting thermal-balance trial point with T_e/T_rad %s: %s",
                candidate[1::2],
                error,
            )
            return np.full_like(candidate, np.nan)

    def calculate_residual_noise(
        self, candidate: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Calculate the Monte Carlo standard error of each residual.

        The residual is evaluated once with each packet batch's estimators;
        the scatter of the batches gives the noise of the residual of the
        summed estimators, including the response of the populations. A
        batch whose ionization balance fails, or whose residuals are not
        finite, is dropped with a warning. Its few packets make that more
        likely than for the summed estimators, and dropping it slightly
        underestimates the noise.

        Parameters
        ----------
        candidate : numpy.ndarray
            Interleaved electron-density fractions and electron-to-radiation
            temperature ratios for each shell.

        Returns
        -------
        numpy.ndarray
            Interleaved standard errors of the electron and heating residuals.

        Raises
        ------
        ValueError
            If fewer than two batches remain.
        """
        batch_residuals = []
        batch_packet_counts = []
        for batch_idx, (evaluator, packet_count) in enumerate(
            zip(self.batch_evaluators(), self.batch_packet_counts, strict=True)
        ):
            try:
                batch_residual = self.calculate_residual(candidate, evaluator)
            except PlasmaIonizationError as error:
                logger.warning(
                    "Dropping packet batch %d from the thermal-balance noise "
                    "estimate: %s",
                    batch_idx,
                    error,
                )
                continue
            if not np.isfinite(batch_residual).all():
                logger.warning(
                    "Dropping packet batch %d from the thermal-balance noise "
                    "estimate: nonfinite residual in shells %s",
                    batch_idx,
                    np.flatnonzero(
                        ~np.isfinite(batch_residual).reshape(-1, 2).all(axis=1)
                    ),
                )
                continue
            batch_residuals.append(batch_residual)
            batch_packet_counts.append(packet_count)
        return calculate_batch_standard_error(
            batch_residuals, batch_packet_counts
        )

    @staticmethod
    def find_converged_shells(
        outer_residuals: npt.NDArray[np.float64],
        outer_residual_limits: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.bool_]:
        """Find the shells whose outer residuals are all within their limits.

        Parameters
        ----------
        outer_residuals : numpy.ndarray
            Interleaved electron and heating residuals for each shell.
        outer_residual_limits : numpy.ndarray
            Limit for each value of ``outer_residuals``.

        Returns
        -------
        numpy.ndarray
            True for each shell with both residuals within their limits.
            Nonfinite residuals are not closed.
        """
        return np.all(
            np.abs(np.asarray(outer_residuals)).reshape(-1, 2)
            <= np.asarray(outer_residual_limits).reshape(-1, 2),
            axis=1,
        )

    @staticmethod
    def build_convergence_stop_options(
        residual_limits: npt.NDArray[np.float64],
    ) -> dict[str, Callable[[OptimizeResult], None]]:
        """Return optimizer options that stop once every shell closes.

        A shell is closed when both of its residuals are within their limits.

        Parameters
        ----------
        residual_limits : numpy.ndarray
            Interleaved electron and heating limits for each shell, in the
            units of the optimizer's residuals.

        Returns
        -------
        dict
            A ``callback`` option for ``scipy.optimize.least_squares``.
        """

        def stop(intermediate_result: OptimizeResult) -> None:
            closed_shells = ElectronThermalBalanceSolver.find_converged_shells(
                intermediate_result.fun, residual_limits
            )
            if closed_shells.all():
                raise StopIteration

        return {"callback": stop}

    def calculate_row_scale(
        self,
        candidate: npt.NDArray[np.float64],
        residual: npt.NDArray[np.float64],
        upper_bound: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Calculate the inverse Jacobian diagonal of each residual row.

        Dividing each row by its own derivative puts the electron and heating
        rows of every shell in candidate units. In inner shells the heating
        row is otherwise ~100 times less sensitive than the electron row, and
        the trust-region solve stalls in their curved valley. Shells are
        independent, so one forward difference per candidate variable, taken
        in every shell at once, gives the whole diagonal.

        Parameters
        ----------
        candidate : numpy.ndarray
            Interleaved electron-density fractions and electron-to-radiation
            temperature ratios for each shell.
        residual : numpy.ndarray
            Interleaved electron and heating residuals at ``candidate``.
        upper_bound : numpy.ndarray
            Upper bounds for each value in ``candidate``. A step that would
            cross a bound is taken backward instead.

        Returns
        -------
        numpy.ndarray
            1 / |d residual_i / d candidate_i| for each residual row, or 1
            where that derivative is zero or not finite.
        """
        row_scale = np.ones_like(candidate)
        for variable_idx in (0, 1):
            values = candidate[variable_idx::2]
            step = ROW_SCALE_RELATIVE_STEP * values
            step = np.where(
                values + step > upper_bound[variable_idx::2], -step, step
            )
            stepped_candidate = candidate.copy()
            stepped_candidate[variable_idx::2] += step
            derivative = (
                self.calculate_trial_residual(stepped_candidate)[
                    variable_idx::2
                ]
                - residual[variable_idx::2]
            ) / step
            usable = np.isfinite(derivative) & (derivative != 0.0)
            if not usable.all():
                logger.warning(
                    "Thermal-balance rows of shells %s keep unit scale: their "
                    "Jacobian diagonal is zero or not finite.",
                    np.flatnonzero(~usable),
                )
            row_scale[variable_idx::2] = np.where(
                usable, 1.0 / np.abs(np.where(usable, derivative, 1.0)), 1.0
            )
        return row_scale

    def solve_multiple_shells(
        self,
        initial_guess: npt.NDArray[np.float64],
        lower_bound: npt.NDArray[np.float64],
        upper_bound: npt.NDArray[np.float64],
        residual_limits: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """
        Solve charge and heating balance for all shells jointly.

        The solve stops once every shell has both residuals within its
        limits, and returns the initial guess unchanged if it already does.
        Each residual row is divided by its Jacobian diagonal at the initial
        guess (see ``calculate_row_scale``). The scaling leaves the roots
        unchanged. The stopping test compares the scaled residuals with the
        limits scaled the same way, which is equivalent to comparing the
        unscaled residuals with the limits.

        Parameters
        ----------
        initial_guess : numpy.ndarray
            Initial interleaved electron-density fractions and temperature
            ratios for each shell.
        lower_bound : numpy.ndarray
            Lower bounds for each value in ``initial_guess``.
        upper_bound : numpy.ndarray
            Upper bounds for each value in ``initial_guess``.
        residual_limits : numpy.ndarray
            Interleaved electron and heating limits for each shell. Zero
            limits never stop the solve early, which gives the optimizer's
            tightest root.

        Returns
        -------
        numpy.ndarray
            Accepted interleaved electron-density fractions and temperature
            ratios.
        """
        initial_residual = self.calculate_trial_residual(initial_guess)
        if self.find_converged_shells(initial_residual, residual_limits).all():
            return initial_guess.copy()
        row_scale = self.calculate_row_scale(
            initial_guess, initial_residual, upper_bound
        )

        def scaled_residual(
            candidate: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            # The optimizer starts at the initial guess, already evaluated.
            if np.array_equal(candidate, initial_guess):
                return initial_residual * row_scale
            return self.calculate_trial_residual(candidate) * row_scale

        no_shells = len(initial_guess) // 2
        return least_squares(
            scaled_residual,
            initial_guess,
            bounds=(lower_bound, upper_bound),
            jac_sparsity=block_diag([np.ones((2, 2))] * no_shells),
            verbose=1,
            method="trf",
            **self.least_squares_settings,
            **self.build_convergence_stop_options(residual_limits * row_scale),
        ).x

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
        Validate the populations and the level and charge closure.

        The electron and heating residuals are not checked: they close at
        their Monte Carlo noise, and open shells are reported by
        ``warn_open_shells``.

        Parameters
        ----------
        evaluation : tardis.plasma.equilibrium.evaluator.PlasmaEquilibriumEvaluation
            Terminal equilibrium state to validate.

        Raises
        ------
        RuntimeError
            If populations or derived states are invalid, or if a level or
            charge residual exceeds its limit in
            ``THERMAL_BALANCE_RESIDUAL_LIMITS``.
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
        Solve all shells jointly and validate the accepted state.

        The closure limits are twice the Monte Carlo standard errors at the
        initial guess. The solve stops once every shell closes; shells still
        open at the evaluation budget are reported by a warning.

        Parameters
        ----------
        initial_guess : numpy.ndarray
            Initial interleaved electron-density fractions and temperature
            ratios for each shell.
        lower_bound : numpy.ndarray
            Lower bounds for each value in ``initial_guess``.
        upper_bound : numpy.ndarray
            Upper bounds for each value in ``initial_guess``.

        Returns
        -------
        accepted_candidate : numpy.ndarray
            Interleaved electron-density fractions and temperature ratios of
            the accepted state.
        evaluation : tardis.plasma.equilibrium.evaluator.PlasmaEquilibriumEvaluation
            Validated terminal equilibrium state for ``accepted_candidate``.

        Raises
        ------
        RuntimeError
            If the populations or the level and charge closure are invalid.
        """
        self.outer_residual_limits = 2.0 * self.calculate_residual_noise(
            initial_guess
        )
        logger.info(
            "Monte Carlo closure limits (2 sigma): electron %s heating %s",
            self.outer_residual_limits[::2],
            self.outer_residual_limits[1::2],
        )
        loose_shells = np.flatnonzero(
            np.any(
                self.outer_residual_limits.reshape(-1, 2)
                > NOISE_LIMITED_CLOSURE_WARNING,
                axis=1,
            )
        )
        if len(loose_shells):
            logger.warning(
                "Monte Carlo noise lets thermal-balance residuals exceed "
                "%.0e in shells %s; consider using more packets.",
                NOISE_LIMITED_CLOSURE_WARNING,
                loose_shells,
            )
        accepted_candidate = self.solve_multiple_shells(
            initial_guess, lower_bound, upper_bound, self.outer_residual_limits
        )
        evaluation = self.evaluate_candidate(accepted_candidate)
        self.validate_evaluation(evaluation)
        self.warn_failed_shells(evaluation)
        return accepted_candidate, evaluation

    def warn_failed_shells(
        self, evaluation: PlasmaEquilibriumEvaluation
    ) -> None:
        """Warn about shells that have not reached their Monte Carlo limits.

        The terminal evaluation is seeded from the converged populations,
        while the optimizer's stopping test uses the frozen level seed. The
        two differ by about the level-closure tolerance, far below the
        Monte Carlo limits.

        Parameters
        ----------
        evaluation : tardis.plasma.equilibrium.evaluator.PlasmaEquilibriumEvaluation
            Terminal equilibrium state.
        """
        outer_residuals = self.calculate_electron_and_thermal_residuals(
            evaluation
        )
        residual_to_limit_ratio = (
            np.abs(outer_residuals) / self.outer_residual_limits
        ).reshape(-1, 2)
        open_shell_idx = np.flatnonzero(
            ~self.find_converged_shells(
                outer_residuals, self.outer_residual_limits
            )
        )
        if len(open_shell_idx):
            logger.warning(
                "Thermal balance did not reach its Monte Carlo closure in "
                "%d of %d shells (|residual| / 2 sigma):\n%s",
                len(open_shell_idx),
                len(residual_to_limit_ratio),
                "\n".join(
                    f"  shell {shell_idx}: electron "
                    f"{residual_to_limit_ratio[shell_idx, 0]:.3g}, heating "
                    f"{residual_to_limit_ratio[shell_idx, 1]:.3g}"
                    for shell_idx in open_shell_idx
                ),
            )
