"""Coordinate normalized electron-density and temperature equilibrium roots."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from time import perf_counter
from typing import Literal, Protocol

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.optimize import OptimizeResult, least_squares
from scipy.sparse import block_diag

logger = logging.getLogger(__name__)


class EquilibriumEvaluator(Protocol):
    """Evaluator interface consumed by :class:`PlasmaEquilibriumSolver`."""

    def evaluate(
        self,
        trial_electron_density: npt.NDArray[np.float64],
        electron_temperature: npt.NDArray[np.float64],
        level_seed: pd.DataFrame,
    ) -> EquilibriumEvaluation:
        """Return electron and fractional-heating residuals."""


class EquilibriumEvaluation(Protocol):
    """Minimum evaluator output required by the normalized root solver."""

    electron_residual: pd.Series | None
    fractional_heating: pd.Series | None
    normalized_population: pd.DataFrame


@dataclass(frozen=True)
class PlasmaEquilibriumSolution:
    """Terminal state, with kinetic closure distinguished from continuation.

    ``maximum_relative_level_residual=None`` denotes unavailable optional
    kinetic diagnostics, not a measured zero. The true default preserves the
    historical minimum evaluator contract for three-argument construction.
    """

    electron_density: npt.NDArray[np.float64]
    electron_temperature: npt.NDArray[np.float64]
    evaluation: EquilibriumEvaluation
    level_converged: bool = True
    maximum_relative_level_residual: float | None = None


class PlasmaEquilibriumSolver:
    """Solve an evaluator in normalized density-temperature coordinates."""

    def __init__(
        self,
        minimum_electron_temperature: float = 1500.0,
        maximum_temperature_ratio: float = 1.5,
        maximum_function_evaluations: int = 100,
        residual_tolerance: float = 1.0e-7,
        fixed_temperature_ratio: float | None = None,
        level_failure_policy: Literal["raise", "warn"] = "raise",
    ) -> None:
        """Configure bounds, solver effort, and optional prescribed temperature.

        Parameters
        ----------
        minimum_electron_temperature : float, optional
            Temperature floor in kelvin for coupled thermal-balance solves.
        maximum_temperature_ratio : float, optional
            Upper T_e/T_rad bound for coupled thermal-balance solves.
        maximum_function_evaluations : int, optional
            Maximum number of nonlinear function evaluations.
        residual_tolerance : float, optional
            Maximum accepted absolute normalized residual.
        fixed_temperature_ratio : float or None, optional
            Prescribe T_e/T_rad and solve only for electron density when set.
            Thermal residuals and temperature bounds are unused in this mode.
            None retains the coupled density-temperature solve.
        level_failure_policy : {"raise", "warn"}, optional
            Reject finite terminal kinetic-closure failure by default, or
            warn and return an explicitly unclosed candidate. Nonfinite or
            physically invalid states and other residual failures stay fatal.
        """
        self.minimum_electron_temperature = minimum_electron_temperature
        self.maximum_temperature_ratio = maximum_temperature_ratio
        self.maximum_function_evaluations = maximum_function_evaluations
        self.residual_tolerance = residual_tolerance
        if fixed_temperature_ratio is not None and (
            not np.isfinite(fixed_temperature_ratio)
            or fixed_temperature_ratio <= 0.0
        ):
            raise ValueError(
                "fixed_temperature_ratio must be finite and positive."
            )
        self.fixed_temperature_ratio = fixed_temperature_ratio
        if level_failure_policy not in ("raise", "warn"):
            raise ValueError("level_failure_policy must be 'raise' or 'warn'.")
        self.level_failure_policy = level_failure_policy

    def solve(
        self,
        evaluator: EquilibriumEvaluator,
        initial_electron_density: npt.ArrayLike,
        initial_electron_temperature: npt.ArrayLike,
        level_seed: pd.DataFrame,
        maximum_electron_density: npt.ArrayLike,
        radiation_temperature: npt.ArrayLike,
    ) -> PlasmaEquilibriumSolution:
        """Solve charge and optionally thermal residuals with fixed rate inputs.

        Intermediate level solves may return physical best-effort candidates.
        If terminal kinetic closure fails, refresh the evaluator's population
        context up to three times and reoptimize charge/temperature as needed.
        All terminal tolerances apply to the exact returned state. The warn
        policy permits only finite kinetic failure after these refinements;
        it never declares the resulting level equilibrium converged.

        Parameters
        ----------
        evaluator : EquilibriumEvaluator
            Evaluator whose rate inputs remain fixed during this solve.
        initial_electron_density : array-like
            Initial electron density in inverse cubic centimeters.
        initial_electron_temperature : array-like
            Initial electron temperature in kelvin.
        level_seed : pandas.DataFrame
            Initial selected-species normalized level populations.
        maximum_electron_density : array-like
            Full-ionization electron-density scale.
        radiation_temperature : array-like
            Radiation-temperature scale for the thermal unknown or the
            prescribed electron temperature.

        Returns
        -------
        PlasmaEquilibriumSolution
            Accepted dimensional state and terminal evaluation.
        """
        initial_density = np.asarray(initial_electron_density, dtype=np.float64)
        initial_temperature = np.asarray(
            initial_electron_temperature, dtype=np.float64
        )
        maximum_density = np.asarray(maximum_electron_density, dtype=np.float64)
        radiation_temperature = np.asarray(
            radiation_temperature, dtype=np.float64
        )
        shell_count = len(initial_density)
        variable_stride = 2 if self.fixed_temperature_ratio is None else 1
        prescribed_temperature = (
            None
            if self.fixed_temperature_ratio is None
            else self.fixed_temperature_ratio * radiation_temperature
        )
        solve_start_time = perf_counter()
        logger.info(
            "Starting plasma equilibrium solve for %d shells (%d variables; "
            "maximum function evaluations %d).",
            shell_count,
            variable_stride * shell_count,
            self.maximum_function_evaluations,
        )
        initial_guess = np.empty(variable_stride * shell_count)
        initial_guess[::variable_stride] = initial_density / maximum_density
        lower_bound = np.empty_like(initial_guess)
        lower_bound[::variable_stride] = np.finfo(float).eps
        upper_bound = np.empty_like(initial_guess)
        upper_bound[::variable_stride] = 1.0
        if prescribed_temperature is None:
            initial_guess[1::2] = initial_temperature / radiation_temperature
            lower_bound[1::2] = (
                self.minimum_electron_temperature / radiation_temperature
            )
            upper_bound[1::2] = self.maximum_temperature_ratio
        initial_guess = np.clip(initial_guess, lower_bound, upper_bound)
        jacobian_sparsity = block_diag(
            [np.ones((variable_stride, variable_stride), dtype=bool)]
            * shell_count,
            format="csr",
        )
        residual_evaluation_count = 0
        optimization_evaluation_start = 0
        previous_optimizer_candidate = initial_guess.copy()
        previous_optimizer_cost = np.nan
        thermal_budgets: dict[bytes, pd.DataFrame] = {}

        def residual(
            candidate: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            nonlocal residual_evaluation_count, previous_optimizer_cost
            residual_evaluation_count += 1
            evaluation_idx = residual_evaluation_count
            evaluation_start_time = perf_counter()
            logger.debug(
                "Starting plasma equilibrium residual evaluation %d.",
                evaluation_idx,
            )
            evaluation = evaluator.evaluate(
                maximum_density * candidate[::variable_stride],
                (
                    radiation_temperature * candidate[1::2]
                    if prescribed_temperature is None
                    else prescribed_temperature
                ),
                level_seed,
            )
            if evaluation.electron_residual is None:
                raise ValueError(
                    "The evaluator must provide electron residuals."
                )
            electron_residual = np.asarray(
                evaluation.electron_residual, dtype=np.float64
            )
            if prescribed_temperature is not None:
                trial_cost = float(
                    0.5 * np.dot(electron_residual, electron_residual)
                )
                if evaluation_idx == optimization_evaluation_start + 1:
                    previous_optimizer_cost = trial_cost
                logger.debug(
                    "Finished plasma equilibrium residual evaluation %d in %.3f s "
                    "(electron-density max %.6e; fixed T_e/T_rad=%.3f).",
                    evaluation_idx,
                    perf_counter() - evaluation_start_time,
                    np.max(np.abs(electron_residual)),
                    self.fixed_temperature_ratio,
                )
                return electron_residual
            if evaluation.fractional_heating is None:
                raise ValueError(
                    "The evaluator must provide electron and thermal residuals."
                )
            result = np.empty_like(candidate)
            thermal_residual = np.asarray(
                evaluation.fractional_heating, dtype=np.float64
            )
            result[::2] = electron_residual
            result[1::2] = thermal_residual
            thermal_budget = getattr(evaluation, "thermal_budget", None)
            if thermal_budget is not None:
                thermal_budgets[candidate.tobytes()] = thermal_budget
            charge_shell_idx = int(np.argmax(np.abs(electron_residual)))
            thermal_shell_idx = int(np.argmax(np.abs(thermal_residual)))
            trial_cost = float(0.5 * np.dot(result, result))
            if evaluation_idx == optimization_evaluation_start + 1:
                previous_optimizer_cost = trial_cost
            logger.debug(
                "Finished plasma equilibrium residual evaluation %d in %.3f "
                "s (electron-density max %.6e in shell %s (signed %.6e, trial n_e "
                "%.6e cm^-3, n_e/n_e,max %.6e); thermal max %.6e in "
                "shell %s (signed %.6e, trial T_e %.6e K, T_e/T_rad %.6e); "
                "trial cost %.6e).",
                evaluation_idx,
                perf_counter() - evaluation_start_time,
                abs(electron_residual[charge_shell_idx]),
                evaluation.electron_residual.index[charge_shell_idx],
                electron_residual[charge_shell_idx],
                maximum_density[charge_shell_idx]
                * candidate[2 * charge_shell_idx],
                candidate[2 * charge_shell_idx],
                abs(thermal_residual[thermal_shell_idx]),
                evaluation.fractional_heating.index[thermal_shell_idx],
                thermal_residual[thermal_shell_idx],
                radiation_temperature[thermal_shell_idx]
                * candidate[2 * thermal_shell_idx + 1],
                candidate[2 * thermal_shell_idx + 1],
                trial_cost,
            )
            return result

        def report_optimizer_step(intermediate_result: OptimizeResult) -> None:
            nonlocal previous_optimizer_candidate, previous_optimizer_cost
            accepted_residual = intermediate_result.fun
            step_norm = float(
                np.max(
                    np.abs(intermediate_result.x - previous_optimizer_candidate)
                )
            )
            cost_reduction = float(
                previous_optimizer_cost - intermediate_result.cost
            )
            charge_shell_idx = int(
                np.argmax(np.abs(accepted_residual[::variable_stride]))
            )
            if prescribed_temperature is not None:
                logger.info(
                    "Outer optimizer iteration %d: cost %.6e, cost reduction "
                    "%.6e, step infinity norm %.6e (nfev=%d); electron max %.3e "
                    "in shell %s (fixed T_e/T_rad=%.3f).",
                    intermediate_result.nit,
                    intermediate_result.cost,
                    cost_reduction,
                    step_norm,
                    intermediate_result.nfev,
                    abs(accepted_residual[charge_shell_idx]),
                    level_seed.columns[charge_shell_idx],
                    self.fixed_temperature_ratio,
                )
                previous_optimizer_candidate = intermediate_result.x.copy()
                previous_optimizer_cost = float(intermediate_result.cost)
                return
            thermal_shell_idx = int(np.argmax(np.abs(accepted_residual[1::2])))
            thermal_shell = level_seed.columns[thermal_shell_idx]
            budget = thermal_budgets.get(intermediate_result.x.tobytes())
            thermal_budgets.clear()
            budget_report = ""
            if budget is not None:
                shell_budget = budget.loc[thermal_shell]
                budget_terms = (
                    ("total_heating", "H"),
                    ("total_cooling", "C"),
                    ("bound_free_heating", "bound-free H"),
                    ("free_bound_cooling", "free-bound C"),
                    ("free_free_heating", "free-free H"),
                    ("free_free_cooling", "free-free C"),
                    (
                        "collisional_ionization_heating",
                        "collisional ionization H",
                    ),
                    (
                        "collisional_ionization_cooling",
                        "collisional ionization C",
                    ),
                    ("collisional_bound_heating", "collisional bound H"),
                    ("collisional_bound_cooling", "collisional bound C"),
                    ("nonthermal_heating", "nonthermal H"),
                )
                budget_report = (
                    "; thermal budget [erg s^-1 cm^-3]: "
                    + ", ".join(
                        f"{label}={shell_budget[name]:.3e}"
                        for name, label in budget_terms
                    )
                )
                heating_deficit = (
                    shell_budget["total_heating"]
                    - shell_budget["total_cooling"]
                )
                budget_report += f", H-C={heating_deficit:.6e}"
                if "electron_density" in shell_budget:
                    budget_report += (
                        f"; thermal state: n_e={shell_budget['electron_density']:.6e}"
                        " cm^-3, "
                        f"T_e={shell_budget['electron_temperature']:.6e} K"
                    )
                if "adiabatic_cooling_estimate" in shell_budget:
                    adiabatic_cooling = shell_budget[
                        "adiabatic_cooling_estimate"
                    ]
                    cooling_ratio = (
                        f"{adiabatic_cooling / heating_deficit:.6e}"
                        if heating_deficit > 0.0
                        else "n/a (H<=C)"
                    )
                    budget_report += (
                        f", t_exp={shell_budget['time_explosion']:.3e} s; "
                        "electron adiabatic estimate (excluded from C)="
                        f"{adiabatic_cooling:.3e} erg s^-1 cm^-3, "
                        f"adiabatic/(H-C)={cooling_ratio}"
                    )
            logger.info(
                "Outer optimizer iteration %d: cost %.6e, cost reduction "
                "%.6e, step infinity norm %.6e (nfev=%d); electron max %.3e "
                "in shell %s; thermal max %.3e in shell %s (signed %.3e, "
                "T_e/T_rad=%.3e)%s.",
                intermediate_result.nit,
                intermediate_result.cost,
                cost_reduction,
                step_norm,
                intermediate_result.nfev,
                abs(accepted_residual[2 * charge_shell_idx]),
                level_seed.columns[charge_shell_idx],
                abs(accepted_residual[2 * thermal_shell_idx + 1]),
                thermal_shell,
                accepted_residual[2 * thermal_shell_idx + 1],
                intermediate_result.x[2 * thermal_shell_idx + 1],
                budget_report,
            )
            previous_optimizer_candidate = intermediate_result.x.copy()
            previous_optimizer_cost = float(intermediate_result.cost)

        def optimize(candidate: npt.NDArray[np.float64]) -> OptimizeResult:
            nonlocal previous_optimizer_candidate, previous_optimizer_cost
            nonlocal optimization_evaluation_start
            optimization_evaluation_start = residual_evaluation_count
            previous_optimizer_candidate = candidate.copy()
            previous_optimizer_cost = np.nan
            thermal_budgets.clear()
            least_squares_start_time = perf_counter()
            result = least_squares(
                residual,
                candidate,
                bounds=(lower_bound, upper_bound),
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
                max_nfev=self.maximum_function_evaluations,
                method="trf",
                jac_sparsity=jacobian_sparsity,
                callback=report_optimizer_step,
            )
            logger.info(
                "Finished plasma equilibrium least-squares optimization in %.3f "
                "s (residual evaluations=%d, nfev=%d, njev=%s, cost=%.3e, "
                "status=%d).",
                perf_counter() - least_squares_start_time,
                residual_evaluation_count,
                result.nfev,
                result.njev,
                result.cost,
                result.status,
            )
            if not result.success:
                raise RuntimeError(
                    "Plasma equilibrium solve failed: " + result.message
                )
            return result

        root = optimize(initial_guess)
        electron_density = maximum_density * root.x[::variable_stride]
        electron_temperature = (
            radiation_temperature * root.x[1::2]
            if prescribed_temperature is None
            else prescribed_temperature
        )
        terminal_validation_start_time = perf_counter()
        seed_evaluation = evaluator.evaluate(
            electron_density,
            electron_temperature,
            level_seed,
        )
        terminal_evaluation = evaluator.evaluate(
            electron_density,
            electron_temperature,
            seed_evaluation.normalized_population,
        )
        terminal_evaluation_count = 2
        level_tolerance = getattr(
            evaluator, "level_residual_tolerance", self.residual_tolerance
        )
        # Three refreshed contexts allow a terminal opacity/ion correction and
        # confirmation while bounding unresolved fixed-point iteration. Each
        # pass retains the existing full inner and outer optimizer effort.
        refine_terminal = getattr(evaluator, "refine_terminal", None)
        refinement_count = 0
        for refinement_pass in range(1, 4):
            relative_levels = getattr(
                terminal_evaluation, "relative_level_residual", None
            )
            if relative_levels is None or refine_terminal is None:
                break
            level_values = np.asarray(relative_levels, dtype=np.float64)
            before_residual = float(np.max(np.abs(level_values)))
            if (
                np.isfinite(level_values).all()
                and before_residual <= level_tolerance
            ):
                break
            level_idx, shell_idx = np.unravel_index(
                np.argmax(np.abs(level_values)), level_values.shape
            )
            logger.info(
                "Starting terminal level refinement %d/3: level %s, shell %s, "
                "relative residual %.3e (tolerance %.3e).",
                refinement_pass,
                relative_levels.index[level_idx],
                relative_levels.columns[shell_idx],
                before_residual,
                level_tolerance,
            )
            refinement_count = refinement_pass
            terminal_evaluation = refine_terminal(
                electron_density, electron_temperature, terminal_evaluation
            )
            terminal_evaluation_count += 1
            level_seed = terminal_evaluation.normalized_population
            outer_residual = np.asarray(
                terminal_evaluation.electron_residual, dtype=np.float64
            )
            if prescribed_temperature is None:
                outer_residual = np.concatenate(
                    [
                        outer_residual,
                        np.asarray(
                            terminal_evaluation.fractional_heating,
                            dtype=np.float64,
                        ),
                    ]
                )
            if not np.isfinite(outer_residual).all() or np.max(
                np.abs(outer_residual)
            ) > min(self.residual_tolerance, 0.1 * level_tolerance):
                # Refreshing physical levels can shift the charge/thermal root.
                # Reoptimize the same coordinates/objective, with no kinetic
                # rows added and the same independent-shell Jacobian pattern.
                candidate = np.empty_like(initial_guess)
                candidate[::variable_stride] = (
                    electron_density / maximum_density
                )
                if prescribed_temperature is None:
                    candidate[1::2] = (
                        electron_temperature / radiation_temperature
                    )
                root = optimize(np.clip(candidate, lower_bound, upper_bound))
                electron_density = maximum_density * root.x[::variable_stride]
                electron_temperature = (
                    radiation_temperature * root.x[1::2]
                    if prescribed_temperature is None
                    else prescribed_temperature
                )
                terminal_evaluation = evaluator.evaluate(
                    electron_density, electron_temperature, level_seed
                )
                terminal_evaluation_count += 1
            after_levels = np.asarray(
                terminal_evaluation.relative_level_residual, dtype=np.float64
            )
            after_residual = float(np.max(np.abs(after_levels)))
            logger.info(
                "Finished terminal level refinement %d/3: relative residual "
                "%.3e -> %.3e; level acceptance=%s.",
                refinement_pass,
                before_residual,
                after_residual,
                np.isfinite(after_levels).all()
                and after_residual <= level_tolerance,
            )
        terminal_residual = np.asarray(
            terminal_evaluation.electron_residual, dtype=np.float64
        )
        if prescribed_temperature is None:
            terminal_residual = np.concatenate(
                [
                    terminal_residual,
                    np.asarray(
                        terminal_evaluation.fractional_heating,
                        dtype=np.float64,
                    ),
                ]
            )
        charge_residual = getattr(terminal_evaluation, "charge_residual", None)
        if charge_residual is not None:
            terminal_residual = np.concatenate(
                [
                    terminal_residual,
                    np.asarray(charge_residual, dtype=np.float64),
                ]
            )
        if not np.all(np.isfinite(terminal_residual)):
            raise RuntimeError(
                "Plasma equilibrium solve returned nonfinite residuals."
            )
        maximum_residual = float(np.max(np.abs(terminal_residual)))
        if maximum_residual > self.residual_tolerance:
            raise RuntimeError(
                "Plasma equilibrium solve did not close its normalized "
                f"residuals: {maximum_residual:.3e} > "
                f"{self.residual_tolerance:.3e}."
            )
        relative_level_residual = getattr(
            terminal_evaluation, "relative_level_residual", None
        )
        # Validate available physical fields without expanding the minimum
        # evaluator protocol. These are arithmetic consistency checks, not
        # looser kinetic acceptance or a population floor.
        population = terminal_evaluation.normalized_population
        for field, values in (
            ("electron_density", electron_density),
            ("electron_temperature", electron_temperature),
            (
                "trial_electron_density",
                getattr(terminal_evaluation, "trial_electron_density", None),
            ),
            (
                "charge_solved_electron_density",
                getattr(
                    terminal_evaluation, "charge_solved_electron_density", None
                ),
            ),
        ):
            if values is None:
                continue
            array = np.asarray(values, dtype=np.float64)
            if not np.isfinite(array).all() or np.any(array <= 0.0):
                raise RuntimeError(
                    f"Plasma equilibrium solve returned invalid {field}."
                )
            if isinstance(values, pd.Series) and not values.index.equals(
                population.columns
            ):
                raise RuntimeError(
                    f"Plasma equilibrium solve returned misaligned {field} shells."
                )
        arithmetic_tolerance = (
            64.0 * np.finfo(float).eps * max(len(population), 1)
        )
        for field in (
            "normalized_population",
            "absolute_level_population",
            "ion_population",
            "tau_sobolev",
            "beta_sobolev",
        ):
            values = getattr(terminal_evaluation, field, None)
            if values is None:
                continue
            array = np.asarray(values, dtype=np.float64)
            if not np.isfinite(array).all() or (
                field != "tau_sobolev" and np.any(array < 0.0)
            ):
                raise RuntimeError(
                    f"Plasma equilibrium solve returned invalid {field}."
                )
            if not values.columns.equals(population.columns):
                raise RuntimeError(
                    f"Plasma equilibrium solve returned misaligned {field} shells."
                )
        if not np.allclose(
            population.sum(axis=0), 1.0, rtol=0.0, atol=arithmetic_tolerance
        ):
            raise RuntimeError(
                "Plasma equilibrium solve returned invalid population normalization."
            )
        ions = getattr(terminal_evaluation, "ion_population", None)
        levels = getattr(terminal_evaluation, "absolute_level_population", None)
        if ions is not None and levels is not None:
            level_groups = levels.groupby(level=["atomic_number", "ion_number"])
            level_totals = level_groups.sum()
            represented_ions = ions.reindex(
                index=level_totals.index, columns=level_totals.columns
            )
            allocation_bound = (level_groups.size() // 2) * np.nextafter(
                0.0, 1.0
            )
            allocation_underflow = level_groups.max().eq(
                0.0
            ) & represented_ions.le(allocation_bound, axis="index")
            consistent = (
                np.isclose(
                    level_totals,
                    represented_ions,
                    rtol=(
                        64.0 * np.finfo(float).eps * level_groups.size()
                    ).to_numpy()[:, None],
                    atol=0.0,
                )
                | allocation_underflow.to_numpy()
            )
            if not consistent.all():
                raise RuntimeError(
                    "Plasma equilibrium solve returned inconsistent level/ion totals."
                )
            species = getattr(evaluator, "equilibrium_species", None)
            if species is not None:
                selected_levels = levels.loc[species]
                selected_population = population
                if isinstance(population.index, pd.MultiIndex):
                    selected_population = population.loc[species]
                expected_levels = selected_population.multiply(
                    ions.loc[species], axis="columns"
                )
                if not selected_levels.index.equals(
                    expected_levels.index
                ) or not np.allclose(
                    selected_levels,
                    expected_levels,
                    rtol=arithmetic_tolerance,
                    atol=0.0,
                ):
                    raise RuntimeError(
                        "Plasma equilibrium solve returned inconsistent selected-level fractions."
                    )
        elemental_density = getattr(evaluator, "elemental_number_density", None)
        if ions is not None and elemental_density is not None:
            nuclei = ions.groupby(level="atomic_number").sum()
            elements = nuclei.index.union(elemental_density.index)
            expected_nuclei = elemental_density.reindex(
                index=elements, columns=nuclei.columns, fill_value=0.0
            )
            nuclei = nuclei.reindex(index=elements, fill_value=0.0)
            stage_counts = (
                ions.groupby(level="atomic_number")
                .size()
                .reindex(elements, fill_value=1)
            )
            if (
                not elemental_density.columns.equals(nuclei.columns)
                or not np.isfinite(np.asarray(expected_nuclei)).all()
                or not np.allclose(
                    nuclei,
                    expected_nuclei,
                    rtol=(64.0 * np.finfo(float).eps * stage_counts).to_numpy()[
                        :, None
                    ],
                    atol=0.0,
                )
            ):
                raise RuntimeError(
                    "Plasma equilibrium solve returned inconsistent conserved nuclei."
                )
        level_converged = True
        maximum_relative_level_residual = None
        if relative_level_residual is not None:
            level_values = np.asarray(relative_level_residual, dtype=np.float64)
            if not np.isfinite(level_values).all():
                raise RuntimeError(
                    "Plasma equilibrium solve returned nonfinite relative level residuals."
                )
            maximum_relative_level_residual = float(
                np.max(np.abs(level_values))
            )
            level_converged = maximum_relative_level_residual <= level_tolerance
            if not level_converged:
                level_idx, shell_idx = np.unravel_index(
                    np.argmax(np.abs(level_values)), level_values.shape
                )
                worst_level = relative_level_residual.index[level_idx]
                if isinstance(worst_level, tuple):
                    worst_level = tuple(
                        value.item() if isinstance(value, np.generic) else value
                        for value in worst_level
                    )
                message = (
                    "Plasma equilibrium solve did not close its relative level "
                    f"residual at level {worst_level}, "
                    f"shell {relative_level_residual.columns[shell_idx]}: "
                    f"{level_values[level_idx, shell_idx]:.3e} "
                    f"(tolerance={level_tolerance:.3e}, "
                    f"terminal refinements={refinement_count})."
                )
                if self.level_failure_policy == "raise":
                    raise RuntimeError(message)
                logger.warning(
                    "%s level_failure_policy=warn; magnitude=%.3e. "
                    "Continuing with this terminal state: NOT converged level equilibrium.",
                    message,
                    maximum_relative_level_residual,
                )
        logger.info(
            "Finished plasma equilibrium terminal validation in %.3f s "
            "after %d total evaluator calls (maximum residual %.3e; total "
            "solve time %.3f s).",
            perf_counter() - terminal_validation_start_time,
            residual_evaluation_count + terminal_evaluation_count,
            maximum_residual,
            perf_counter() - solve_start_time,
        )
        return PlasmaEquilibriumSolution(
            electron_density,
            electron_temperature,
            terminal_evaluation,
            level_converged,
            maximum_relative_level_residual,
        )
