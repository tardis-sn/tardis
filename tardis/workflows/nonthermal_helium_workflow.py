"""Compose high-energy transport with coupled nonthermal helium equilibrium."""

from __future__ import annotations

import logging
from collections.abc import Callable, Collection, Mapping
from dataclasses import dataclass
from typing import Protocol

import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.energy_input.gamma_ray_transport import (
    calculate_ejecta_velocity_volume,
)
from tardis.io.atom_data.parse_atom_data import parse_atom_data
from tardis.opacities.continuum.continuum_state import ContinuumOpacityState
from tardis.opacities.continuum.macro_atom_state import ContinuumMacroAtomState
from tardis.opacities.macro_atom.macroatom_solver import (
    ContinuumMacroAtomSolver,
)
from tardis.opacities.opacity_solver import OpacitySolver
from tardis.opacities.tau_sobolev import (
    SOBOLEV_COEFFICIENT,
    calculate_beta_sobolev,
    calculate_sobolev_line_opacity,
)
from tardis.plasma.assembly import PlasmaSolverFactory
from tardis.plasma.equilibrium.evaluator import (
    PlasmaEquilibriumEvaluator,
    calculate_lte_populations,
)
from tardis.plasma.equilibrium.inputs import (
    ContinuumCoefficientState,
    DiluteRadiationInputs,
    FixedRateState,
    ShellNumberDensity,
    SobolevInputs,
)
from tardis.plasma.equilibrium.ion_populations import IonPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import (
    AnalyticIonRateMatrix,
    EstimatedIonRateMatrix,
    RateMatrix,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    EstimatedPhotoionizationRateSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.equilibrium.rates.util import (
    reindex_ion_population_to_level_population,
)
from tardis.plasma.equilibrium.solver import (
    PlasmaEquilibriumSolution,
    PlasmaEquilibriumSolver,
)
from tardis.plasma.equilibrium.thermal_balance import ThermalBalanceSolver
from tardis.plasma.nonthermal.pynonthermal_solver import (
    NonthermalPlasmaInput,
    NonthermalPlasmaState,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.transport.montecarlo.configuration import montecarlo_globals
from tardis.transport.montecarlo.estimators import init_estimators_continuum
from tardis.transport.montecarlo.modes.iip.solver import (
    MCTransportSolverIIP,
)
from tardis.workflows.high_energy.deposition import (
    DepositionState,
    prepare_deposition_state,
)
from tardis.workflows.high_energy.tardis_he_workflow import (
    TARDISHEWorkflow,
    TARDISHEWorkflowResult,
)
from tardis.workflows.nonthermal_helium_diagnostics import (
    canonicalize_underflowed_populations,
)
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow

logger = logging.getLogger(__name__)


class NonthermalRateSolver(Protocol):
    """Nonthermal rate provider used by the coupled plasma solver."""

    def solve(
        self, plasma_input: NonthermalPlasmaInput
    ) -> NonthermalPlasmaState:
        """Return fixed rates for one candidate plasma state."""


class OpticalPlasma(Protocol):
    """Mutable optical-plasma fields required by helium coupling."""

    electron_densities: pd.Series
    t_electrons: npt.ArrayLike
    ion_number_density: pd.DataFrame
    level_number_density: pd.DataFrame
    level_boltzmann_factor: pd.DataFrame
    partition_function: pd.DataFrame
    number_density: pd.DataFrame

    def update(self, **kwargs: object) -> None:
        """Commit an accepted plasma state."""


class OpticalSimulationState(Protocol):
    """Optical model fields required by helium coupling."""

    t_radiative: npt.ArrayLike
    time_explosion: u.Quantity
    v_inner: u.Quantity
    v_outer: u.Quantity


class OpticalWorkflow(Protocol):
    """Optical workflow state consumed by the coupled plasma solver."""

    plasma_solver: OpticalPlasma
    simulation_state: OpticalSimulationState

    def run(self) -> None:
        """Run the optical Monte Carlo iteration."""

    def build_helium_evaluator(
        self,
        plasma_state: CoupledPlasmaState,
        fixed_rates: FixedRateState,
        bound_bound_species_filter: tuple[int, int] | None,
    ) -> PlasmaEquilibriumEvaluator:
        """Build an evaluator from the current continuum estimators."""


@dataclass(frozen=True)
class HighEnergySettings:
    """Arguments for one frozen high-energy transport calculation."""

    time_start: float
    time_end: float
    number_of_packets: int
    time_steps: int
    time_space: str
    seed: int
    positronium_fraction: float
    spectrum_bins: int
    grey_opacity: float = -1.0

    def as_run_arguments(self) -> dict[str, object]:
        """Return names accepted by :meth:`TARDISHEWorkflow.run`."""
        return {
            "time_start": self.time_start,
            "time_end": self.time_end,
            "number_of_packets": self.number_of_packets,
            "time_steps": self.time_steps,
            "time_space": self.time_space,
            "seed": self.seed,
            "fp": self.positronium_fraction,
            "spectrum_bins": self.spectrum_bins,
            "grey_opacity": self.grey_opacity,
        }


@dataclass(frozen=True)
class NonthermalCouplingSettings:
    """Convergence and relaxation settings for the outer Picard iteration.

    Level populations satisfy ``abs(new - old) <= atol + rtol *
    max(abs(new), abs(old))`` componentwise. Here ``atol`` is
    ``absolute_population_tolerance`` in cm^-3, not a denominator floor.
    Other quantities retain their relative or fraction-change criteria.
    ``relaxation`` damps electron density and temperature and initializes ion
    damping. Adaptive ion damping is enabled by default, with independent
    factors per element/shell, bounds [0.25, 1], and two stable comparisons at
    an elemental-fraction tolerance of 1e-3 before each promotion. New
    within-ion level fractions are applied without damping and normalized to
    the damped ion totals.
    """

    maximum_iterations: int = 20
    relative_tolerance: float = 1.0e-3
    absolute_population_tolerance: float = 1.0e-20
    relaxation: float = 0.5
    adaptive_ion_relaxation: bool = True
    minimum_ion_relaxation: float = 0.25
    maximum_ion_relaxation: float = 1.0
    ion_candidate_stability_tolerance: float = 1.0e-3
    ion_stable_comparisons: int = 2

    def __post_init__(self) -> None:
        """Validate the bounded adaptive-ion controller's settings."""
        if not self.adaptive_ion_relaxation:
            return
        if not (
            0.0
            < self.minimum_ion_relaxation
            <= self.relaxation
            <= self.maximum_ion_relaxation
            <= 1.0
        ):
            raise ValueError(
                "Adaptive ion relaxation requires 0 < minimum <= relaxation "
                "<= maximum <= 1."
            )
        if (
            not np.isfinite(self.ion_candidate_stability_tolerance)
            or self.ion_candidate_stability_tolerance <= 0.0
            or self.ion_stable_comparisons < 1
        ):
            raise ValueError(
                "Adaptive ion relaxation requires a positive finite candidate "
                "stability tolerance and at least one stable comparison."
            )


class _IonRelaxationController:
    """Adapt one convex mixing factor per element/shell from candidate history."""

    def __init__(
        self,
        settings: NonthermalCouplingSettings,
        elemental_number_density: pd.DataFrame,
    ) -> None:
        """Initialize candidate history independently for each coupling solve."""
        self.settings = settings
        self.elemental_number_density = elemental_number_density.copy(deep=True)
        self.factors = pd.DataFrame(
            settings.relaxation,
            index=elemental_number_density.index,
            columns=elemental_number_density.columns,
        )
        self.stable_counts = (self.factors * 0).astype(int)
        self.previous_fractions: pd.DataFrame | None = None
        self.previous_increment: pd.DataFrame | None = None
        self.previous_change: pd.DataFrame | None = None
        self.diagnostics = pd.DataFrame()

    def update(self, candidate: pd.DataFrame) -> pd.DataFrame:
        """Return bounded factors without changing candidate densities or totals.

        Stability uses maximum candidate-to-candidate changes in elemental
        number fraction. Material reversals or more-than-twofold growth
        halve the factor; two stable comparisons promote it by 0.25. Small
        direction reversals below the stability threshold are ignored.
        """
        settings = self.settings
        reasons = pd.DataFrame(
            "fixed", index=self.factors.index, columns=self.factors.columns
        )
        change = self.factors * np.nan
        if settings.adaptive_ion_relaxation:
            atomic_numbers = candidate.index.get_level_values("atomic_number")
            denominators = self.elemental_number_density.reindex(
                index=atomic_numbers, columns=candidate.columns
            ).to_numpy(dtype=float)
            values = candidate.to_numpy(dtype=float)
            if (
                not np.isfinite(values).all()
                or not np.isfinite(denominators).all()
                or np.any(values < 0.0)
                or np.any((denominators == 0.0) & (values != 0.0))
            ):
                raise ValueError(
                    "Adaptive ion relaxation requires finite nonnegative ion "
                    "densities with represented, nonzero-abundance parents."
                )
            fractions = pd.DataFrame(
                np.divide(
                    values,
                    denominators,
                    out=np.zeros_like(values),
                    where=denominators > 0.0,
                ),
                index=candidate.index,
                columns=candidate.columns,
            )
            reasons[:] = "initial"
            if self.previous_fractions is not None:
                increment = fractions - self.previous_fractions
                change = (
                    increment.abs()
                    .groupby(level="atomic_number")
                    .max()
                    .reindex(
                        index=self.factors.index,
                        columns=self.factors.columns,
                        fill_value=0.0,
                    )
                )
                stable = change <= settings.ion_candidate_stability_tolerance
                self.stable_counts = (self.stable_counts + 1).where(stable, 0)
                self.stable_counts = self.stable_counts.clip(
                    upper=settings.ion_stable_comparisons
                )
                backoff = pd.DataFrame(
                    False,
                    index=self.factors.index,
                    columns=self.factors.columns,
                )
                if (
                    self.previous_increment is not None
                    and self.previous_change is not None
                ):
                    dot_product = (
                        (increment * self.previous_increment)
                        .groupby(level="atomic_number")
                        .sum()
                        .reindex(
                            index=self.factors.index,
                            columns=self.factors.columns,
                            fill_value=0.0,
                        )
                    )
                    material = (
                        change > settings.ion_candidate_stability_tolerance
                    )
                    backoff = material & (
                        (change > 2.0 * self.previous_change)
                        | (
                            (dot_product < 0.0)
                            & (
                                self.previous_change
                                > settings.ion_candidate_stability_tolerance
                            )
                        )
                    )
                promote = (
                    (self.stable_counts >= settings.ion_stable_comparisons)
                    & (self.factors < settings.maximum_ion_relaxation)
                    & ~backoff
                )
                self.factors = self.factors.mask(
                    promote,
                    (self.factors + 0.25).clip(
                        upper=settings.maximum_ion_relaxation
                    ),
                )
                self.factors = self.factors.mask(
                    backoff,
                    (0.5 * self.factors).clip(
                        lower=settings.minimum_ion_relaxation
                    ),
                )
                self.stable_counts = self.stable_counts.mask(
                    promote | backoff, 0
                )
                reasons[:] = "hold"
                reasons = (
                    reasons.mask(stable, "stable")
                    .mask(promote, "promote")
                    .mask(backoff, "backoff")
                )
                self.previous_increment = increment
                self.previous_change = change
            self.previous_fractions = fractions
            absent = self.elemental_number_density == 0.0
            self.factors = self.factors.mask(absent, 1.0)
            self.stable_counts = self.stable_counts.mask(absent, 0)
            reasons = reasons.mask(absent, "zero_abundance")
        self.diagnostics = pd.DataFrame(
            {
                "factor": self.factors.to_numpy().ravel(),
                "stable_comparisons": self.stable_counts.to_numpy().ravel(),
                "candidate_change": change.to_numpy().ravel(),
                "reason": reasons.to_numpy().ravel(),
            },
            index=pd.MultiIndex.from_product(
                [self.factors.index, self.factors.columns],
                names=["atomic_number", "shell"],
            ),
        )
        return self.factors.copy(deep=True)


@dataclass(frozen=True)
class CoupledPlasmaState:
    """Candidate plasma fields refreshed between nonthermal solves."""

    electron_density: npt.NDArray[np.float64]
    electron_temperature: npt.NDArray[np.float64]
    ion_number_density: pd.DataFrame
    level_number_density: pd.DataFrame


@dataclass(frozen=True)
class NonthermalCouplingResult:
    """Converged plasma and retained raw equilibrium diagnostics.

    Attributes
    ----------
    plasma_state : CoupledPlasmaState
        Accepted coupling state, including allocation-underflow corrections.
        These same populations are committed to the optical plasma.
    nonthermal_state : NonthermalPlasmaState
        Frozen rates used by the terminal equilibrium solve.
    equilibrium_solution : PlasmaEquilibriumSolution
        Raw terminal evaluator result, retained without mutation. It may
        contain a subnormal parent with all-zero levels that was canonicalized
        to zero in `plasma_state`.
    convergence_history : pandas.DataFrame
        Component changes for each nonthermal coupling iteration.
    """

    plasma_state: CoupledPlasmaState
    nonthermal_state: NonthermalPlasmaState
    equilibrium_solution: PlasmaEquilibriumSolution
    convergence_history: pd.DataFrame


EvaluatorFactory = Callable[
    [
        CoupledPlasmaState,
        FixedRateState,
        OpticalWorkflow,
        tuple[int, int] | None,
    ],
    PlasmaEquilibriumEvaluator,
]


class NonthermalHeliumPlasmaSolver:
    """Iterate pynonthermal outside frozen evaluator equilibrium solves.

    Attributes
    ----------
    last_convergence_history : pandas.DataFrame
        Component changes from the latest coupling attempt, including failure.
    last_iteration_diagnostics : pandas.DataFrame
        Worst entry for each compared quantity, measured before relaxation.
    last_zero_target_levels : pandas.DataFrame
        Failing level/shell entries whose nonnegative candidate population is
        at or below the absolute tolerance. Includes the predicted relaxed
        population and its ratio to the current population, not a new solve.
    last_input_state : CoupledPlasmaState or None
        Input to the latest completed equilibrium solve.
    last_candidate_state : CoupledPlasmaState or None
        Unrelaxed output of the latest completed equilibrium solve, after
        allocation-underflow canonicalization. Raw evaluator tables are not
        mutated. A failed coupling attempt does not commit this candidate to
        the optical plasma.
    last_ion_relaxation_history : pandas.DataFrame
        Factors, counters, candidate changes, and adjustment reasons indexed
        by iteration, atomic number, and shell for the latest coupling attempt.
    """

    def __init__(
        self,
        nonthermal_solver: NonthermalRateSolver,
        evaluator_factory: EvaluatorFactory | None = None,
        equilibrium_solver: PlasmaEquilibriumSolver | None = None,
        settings: NonthermalCouplingSettings = NonthermalCouplingSettings(),
        bound_bound_species_filter: tuple[int, int] | None = (2, 0),
    ) -> None:
        """Configure coupled nonthermal and equilibrium owners.

        Parameters
        ----------
        nonthermal_solver : NonthermalRateSolver
            Provider of fixed nonthermal rates for candidate plasma states.
        evaluator_factory : EvaluatorFactory, optional
            Override builder receiving the candidate state, fixed rates,
            optical workflow, and optional bound-bound species filter. When
            omitted, the optical workflow owns evaluator construction.
        equilibrium_solver : PlasmaEquilibriumSolver, optional
            Normalized charge solver. The default prescribes T_e = 0.9 T_rad.
            Pass PlasmaEquilibriumSolver() to restore coupled thermal balance.
        settings : NonthermalCouplingSettings, optional
            Picard convergence and relaxation settings.
        bound_bound_species_filter : tuple[int, int] or None, optional
            Atomic number and zero-based ion number whose bound-bound inputs
            the evaluator factory should retain. ``None`` requests no line
            filtering. The helium workflow defaults to He I, ``(2, 0)``.
        """
        self.nonthermal_solver = nonthermal_solver
        self.evaluator_factory = evaluator_factory
        self.equilibrium_solver = equilibrium_solver or PlasmaEquilibriumSolver(
            fixed_temperature_ratio=0.9
        )
        self.bound_bound_species_filter = bound_bound_species_filter
        self.settings = settings
        self.deposition_state: DepositionState | None = None
        self.last_result: NonthermalCouplingResult | None = None
        self.last_convergence_history = pd.DataFrame()
        self.last_iteration_diagnostics = pd.DataFrame()
        self.last_zero_target_levels = pd.DataFrame()
        self.last_ion_relaxation_history = pd.DataFrame()
        self.last_input_state: CoupledPlasmaState | None = None
        self.last_candidate_state: CoupledPlasmaState | None = None

    @staticmethod
    def _maximum_relative_change(
        candidate: npt.ArrayLike,
        previous: npt.ArrayLike,
        floor: float,
    ) -> float:
        """Return a floored componentwise relative change."""
        candidate = np.asarray(candidate, dtype=np.float64)
        previous = np.asarray(previous, dtype=np.float64)
        return float(
            np.max(
                np.abs(candidate - previous)
                / np.maximum(np.abs(previous), floor)
            )
        )

    @staticmethod
    def _calculate_change_diagnostic(
        candidate: pd.Series | pd.DataFrame,
        current: pd.Series | pd.DataFrame,
        floor: float | None,
        unit: str,
        change_values: npt.NDArray[np.float64] | None = None,
    ) -> dict[str, object]:
        """Locate the largest pre-relaxation change using the convergence metric.

        ``change_values`` supplies the mixed-tolerance level metric; otherwise
        the original absolute or floored relative comparison is used.
        """
        candidate_values = candidate.to_numpy(dtype=float)
        current_values = current.to_numpy(dtype=float)
        absolute_change = np.abs(candidate_values - current_values)
        if change_values is not None:
            change = change_values
        elif floor is None:
            change = absolute_change
        else:
            change = absolute_change / np.maximum(np.abs(current_values), floor)
        position = np.unravel_index(np.argmax(change), change.shape)
        if isinstance(candidate, pd.DataFrame):
            row_idx, shell_idx = position
            state_index = candidate.index[row_idx]
            if isinstance(state_index, tuple):
                state_index = tuple(
                    value.item() if isinstance(value, np.generic) else value
                    for value in state_index
                )
            elif isinstance(state_index, np.generic):
                state_index = state_index.item()
            shell = candidate.columns[shell_idx]
        else:
            state_index = None
            shell = candidate.index[position[0]]
        if isinstance(shell, np.generic):
            shell = shell.item()
        return {
            "change": float(change[position]),
            "shell": shell,
            "state_index": state_index,
            "current": float(current_values[position]),
            "candidate": float(candidate_values[position]),
            "absolute_change": float(absolute_change[position]),
            "floor": floor if change_values is None else None,
            "absolute_tolerance": floor if change_values is not None else None,
            "unit": unit,
        }

    @staticmethod
    def _rate_change(
        candidate: FixedRateState,
        previous: FixedRateState | None,
        include_heating: bool = True,
    ) -> float:
        """Return the largest ionization change, optionally including heating."""
        if previous is None:
            return np.inf
        changes = []
        for candidate_rate, previous_rate in (
            (candidate.ionization_rates, previous.ionization_rates),
            (
                candidate.heating_rate_density if include_heating else None,
                previous.heating_rate_density if include_heating else None,
            ),
        ):
            if candidate_rate is None and previous_rate is None:
                continue
            if candidate_rate is None or previous_rate is None:
                return np.inf
            changes.append(
                NonthermalHeliumPlasmaSolver._maximum_relative_change(
                    candidate_rate,
                    previous_rate,
                    np.finfo(float).tiny,
                )
            )
        return max(changes, default=0.0)

    @staticmethod
    def _check_finite_level_populations(
        levels: pd.DataFrame, stage: str
    ) -> None:
        """Stop at the first nonfinite level, identifying its update boundary."""
        bad = np.argwhere(~np.isfinite(levels.to_numpy(dtype=float)))
        if len(bad):
            row_idx, shell_idx = bad[0]
            level = tuple(int(value) for value in levels.index[row_idx])
            shell = levels.columns[shell_idx]
            raise ValueError(
                f"Nonfinite {stage} level population: level={level}, "
                f"shell={shell}, value={levels.iloc[row_idx, shell_idx]}."
            )

    def _initial_state(
        self, optical_workflow: OpticalWorkflow
    ) -> CoupledPlasmaState:
        """Copy the current optical plasma state."""
        plasma = optical_workflow.plasma_solver
        return CoupledPlasmaState(
            np.asarray(plasma.electron_densities, dtype=np.float64),
            np.asarray(plasma.t_electrons, dtype=np.float64),
            plasma.ion_number_density.copy(deep=True),
            plasma.level_number_density.copy(deep=True),
        )

    def solve(
        self, optical_workflow: OpticalWorkflow
    ) -> NonthermalCouplingResult:
        """Converge and commit one optical iteration's helium plasma state."""
        if self.deposition_state is None:
            raise RuntimeError("High-energy deposition has not been prepared.")
        self.last_convergence_history = pd.DataFrame()
        self.last_iteration_diagnostics = pd.DataFrame()
        self.last_zero_target_levels = pd.DataFrame()
        self.last_ion_relaxation_history = pd.DataFrame()
        self.last_input_state = None
        self.last_candidate_state = None
        current = self._initial_state(optical_workflow)
        self.last_input_state = current
        self._check_finite_level_populations(
            current.level_number_density, "initial"
        )
        initial_levels, initial_ions = canonicalize_underflowed_populations(
            current.level_number_density, current.ion_number_density
        )
        current = CoupledPlasmaState(
            current.electron_density,
            current.electron_temperature,
            initial_ions,
            initial_levels,
        )
        plasma = optical_workflow.plasma_solver
        shell_columns = plasma.number_density.columns
        logger.info(
            "Solving coupled nonthermal helium plasma in all %d shells: %s",
            len(shell_columns),
            shell_columns.tolist(),
        )
        logger.info(
            "Selective relaxation: level fractions are undamped; "
            "electron density and temperature use relaxation=%.3f; "
            "adaptive ion relaxation=%s.",
            self.settings.relaxation,
            self.settings.adaptive_ion_relaxation,
        )
        maximum_electron_density = (
            plasma.number_density.multiply(
                plasma.number_density.index.to_numpy(dtype=float), axis=0
            )
            .sum()
            .to_numpy(dtype=float)
        )
        complete_radiation_temperature = np.asarray(
            optical_workflow.simulation_state.t_radiative, dtype=np.float64
        )
        temperature_ratio = self.equilibrium_solver.fixed_temperature_ratio
        if temperature_ratio is not None:
            current = CoupledPlasmaState(
                current.electron_density,
                temperature_ratio * complete_radiation_temperature,
                current.ion_number_density,
                current.level_number_density,
            )
            logger.info(
                "Prescribing electron temperature T_e/T_rad=%.3f; "
                "thermal balance is disabled for this helium solve.",
                temperature_ratio,
            )
        previous_fixed_rates = None
        convergence_rows = []
        ion_controller = _IonRelaxationController(
            self.settings, plasma.number_density.loc[:, shell_columns]
        )
        ion_relaxation_rows = []
        terminal_nonthermal = None
        terminal_solution = None

        for iteration in range(1, self.settings.maximum_iterations + 1):
            nonthermal_state = self.nonthermal_solver.solve(
                NonthermalPlasmaInput(
                    self.deposition_state.rate_density,
                    current.electron_density / u.cm**3,
                    current.electron_temperature * u.K,
                    current.ion_number_density,
                )
            )
            if self.evaluator_factory is None:
                evaluator = optical_workflow.build_helium_evaluator(
                    current,
                    nonthermal_state.fixed_rates,
                    self.bound_bound_species_filter,
                )
            else:
                evaluator = self.evaluator_factory(
                    current,
                    nonthermal_state.fixed_rates,
                    optical_workflow,
                    self.bound_bound_species_filter,
                )
            if temperature_ratio is not None:
                evaluator.thermal_balance_solver = None
            he_i_level_seed = current.level_number_density.loc[(2, 0)].copy(
                deep=True
            )
            level_seed = he_i_level_seed.divide(
                he_i_level_seed.sum(axis=0), axis="columns"
            )
            if not np.isfinite(level_seed.to_numpy(dtype=float)).all():
                raise ValueError(
                    "He I level populations cannot provide a normalized "
                    "equilibrium seed."
                )
            solution = self.equilibrium_solver.solve(
                evaluator,
                current.electron_density,
                current.electron_temperature,
                level_seed,
                maximum_electron_density,
                complete_radiation_temperature,
            )
            evaluation = solution.evaluation
            if evaluation.ion_population is None:
                raise RuntimeError(
                    "The evaluator did not solve ion populations."
                )
            candidate = CoupledPlasmaState(
                solution.electron_density,
                solution.electron_temperature,
                evaluation.ion_population,
                evaluation.absolute_level_population,
            )
            self.last_input_state = current
            self.last_candidate_state = candidate
            self._check_finite_level_populations(
                candidate.level_number_density, "candidate"
            )
            candidate_levels, candidate_ions = (
                canonicalize_underflowed_populations(
                    candidate.level_number_density, candidate.ion_number_density
                )
            )
            candidate = CoupledPlasmaState(
                candidate.electron_density,
                candidate.electron_temperature,
                candidate_ions,
                candidate_levels,
            )
            self.last_candidate_state = candidate
            relaxation = self.settings.relaxation
            ion_factors = ion_controller.update(candidate.ion_number_density)
            ion_relaxation_rows.append(
                ion_controller.diagnostics.copy(deep=True)
            )
            self.last_ion_relaxation_history = pd.concat(
                ion_relaxation_rows,
                keys=range(1, iteration + 1),
                names=["iteration"],
            )
            stage_factors = pd.DataFrame(
                ion_factors.reindex(
                    index=candidate.ion_number_density.index.get_level_values(
                        "atomic_number"
                    ),
                    columns=shell_columns,
                ).to_numpy(),
                index=candidate.ion_number_density.index,
                columns=shell_columns,
            )
            relaxed_ions = (
                stage_factors * candidate.ion_number_density
                + (1.0 - stage_factors) * current.ion_number_density
            )
            candidate_level_totals = candidate.level_number_density.groupby(
                level=["atomic_number", "ion_number"]
            ).transform("sum")
            current_level_totals = current.level_number_density.groupby(
                level=["atomic_number", "ion_number"]
            ).transform("sum")
            candidate_fractions = candidate.level_number_density.div(
                candidate_level_totals.where(candidate_level_totals != 0.0)
            )
            # A vanished candidate ion supplies no new level distribution.
            # Keep its previous fractions while its damped population survives.
            candidate_fractions = candidate_fractions.mask(
                candidate_level_totals == 0.0,
                current.level_number_density.div(
                    current_level_totals.where(current_level_totals != 0.0)
                ),
            )
            relaxed_ion_totals = reindex_ion_population_to_level_population(
                relaxed_ions, candidate.level_number_density, next_higher=False
            )
            relaxed_levels = (candidate_fractions * relaxed_ion_totals).mask(
                relaxed_ion_totals == 0.0, 0.0
            )
            relaxed_levels, relaxed_ions = canonicalize_underflowed_populations(
                relaxed_levels, relaxed_ions
            )
            self._check_finite_level_populations(relaxed_levels, "relaxed")
            relaxed_state = CoupledPlasmaState(
                relaxation * candidate.electron_density
                + (1.0 - relaxation) * current.electron_density,
                relaxation * candidate.electron_temperature
                + (1.0 - relaxation) * current.electron_temperature,
                relaxed_ions,
                relaxed_levels,
            )
            density_change = self._maximum_relative_change(
                candidate.electron_density,
                current.electron_density,
                np.finfo(float).tiny,
            )
            temperature_change = self._maximum_relative_change(
                candidate.electron_temperature,
                current.electron_temperature,
                np.finfo(float).tiny,
            )
            candidate_helium = candidate.ion_number_density.xs(
                2, level="atomic_number"
            )
            current_helium = current.ion_number_density.xs(
                2, level="atomic_number"
            )
            candidate_he_ii_fraction = candidate_helium.loc[1] / (
                candidate_helium.sum(axis=0)
            )
            current_he_ii_fraction = current_helium.loc[1] / (
                current_helium.sum(axis=0)
            )
            he_ii_fraction_change = float(
                np.max(
                    np.abs(candidate_he_ii_fraction - current_he_ii_fraction)
                )
            )
            nonthermal_change = self._rate_change(
                nonthermal_state.fixed_rates,
                previous_fixed_rates,
                include_heating=temperature_ratio is None,
            )
            other_maximum_change = max(
                density_change,
                temperature_change,
                he_ii_fraction_change,
                nonthermal_change,
            )
            current_levels = current.level_number_density.to_numpy(dtype=float)
            candidate_levels = candidate.level_number_density.to_numpy(
                dtype=float
            )
            level_absolute_change = np.abs(candidate_levels - current_levels)
            level_allowed_change = (
                self.settings.absolute_population_tolerance
                + self.settings.relative_tolerance
                * np.maximum(np.abs(current_levels), np.abs(candidate_levels))
            )
            # Keep the reported score on the common rtol threshold. In the
            # populated limit this tends to the symmetric relative change.
            level_scaled_change = self.settings.relative_tolerance * np.divide(
                level_absolute_change,
                level_allowed_change,
                out=np.zeros_like(level_absolute_change),
                where=level_allowed_change != 0.0,
            )
            level_change = float(np.max(level_scaled_change))
            maximum_change = max(other_maximum_change, level_change)
            failing_levels = (level_absolute_change > level_allowed_change) | (
                ~np.isfinite(level_scaled_change)
            )
            zero_target_mask = (
                failing_levels
                & (candidate_levels >= 0.0)
                & (
                    candidate_levels
                    <= self.settings.absolute_population_tolerance
                )
                & (current_levels > 0.0)
            )
            row_positions, shell_positions = np.nonzero(zero_target_mask)
            target_current = current_levels[zero_target_mask]
            target_candidate = candidate_levels[zero_target_mask]
            relaxed_population = relaxed_levels.to_numpy(dtype=float)[
                zero_target_mask
            ]
            self.last_zero_target_levels = pd.DataFrame(
                {
                    "state_index": list(
                        candidate.level_number_density.index[row_positions]
                    ),
                    "shell": candidate.level_number_density.columns[
                        shell_positions
                    ],
                    "current": target_current,
                    "candidate": target_candidate,
                    "absolute_change": level_absolute_change[zero_target_mask],
                    "scaled_change": level_scaled_change[zero_target_mask],
                    "allowed_change": level_allowed_change[zero_target_mask],
                    "relaxed_population": relaxed_population,
                    "relaxed_to_current_ratio": (
                        relaxed_population / target_current
                    ),
                }
            ).sort_values("absolute_change", ascending=False)
            zero_target_count = len(self.last_zero_target_levels)
            zero_target_absolute_change = float(
                np.max(level_absolute_change[zero_target_mask], initial=0.0)
            )
            convergence_rows.append(
                {
                    "iteration": iteration,
                    "electron_density": density_change,
                    "electron_temperature": temperature_change,
                    "he_ii_fraction": he_ii_fraction_change,
                    "level_population": level_change,
                    "nonthermal_rates": nonthermal_change,
                    "maximum_change": maximum_change,
                    "failing_level_count": int(
                        np.count_nonzero(failing_levels)
                    ),
                    "zero_target_level_count": zero_target_count,
                    "zero_target_maximum_absolute_change": (
                        zero_target_absolute_change
                    ),
                }
            )
            self.last_convergence_history = pd.DataFrame(
                convergence_rows
            ).set_index("iteration")
            comparisons = {
                "electron_density": (
                    pd.Series(candidate.electron_density, index=shell_columns),
                    pd.Series(current.electron_density, index=shell_columns),
                    np.finfo(float).tiny,
                    "cm^-3",
                ),
                "electron_temperature": (
                    pd.Series(
                        candidate.electron_temperature, index=shell_columns
                    ),
                    pd.Series(
                        current.electron_temperature, index=shell_columns
                    ),
                    np.finfo(float).tiny,
                    "K",
                ),
                "he_ii_fraction": (
                    candidate_he_ii_fraction,
                    current_he_ii_fraction,
                    None,
                    "dimensionless",
                ),
                "level_population": (
                    candidate.level_number_density,
                    current.level_number_density,
                    self.settings.absolute_population_tolerance,
                    "cm^-3",
                    level_scaled_change,
                ),
            }
            if previous_fixed_rates is not None:
                for name, new_rate, old_rate, unit in (
                    (
                        "nonthermal_ionization",
                        nonthermal_state.fixed_rates.ionization_rates,
                        previous_fixed_rates.ionization_rates,
                        "s^-1",
                    ),
                    (
                        "nonthermal_heating",
                        nonthermal_state.fixed_rates.heating_rate_density,
                        previous_fixed_rates.heating_rate_density,
                        "erg s^-1 cm^-3",
                    ),
                ):
                    if (
                        name == "nonthermal_heating"
                        and temperature_ratio is not None
                    ):
                        continue
                    if (
                        new_rate is not None
                        and old_rate is not None
                        and new_rate.size
                    ):
                        comparisons[name] = (
                            new_rate,
                            old_rate,
                            np.finfo(float).tiny,
                            unit,
                        )
            diagnostics = [
                {
                    "quantity": name,
                    **self._calculate_change_diagnostic(*arguments),
                }
                for name, arguments in comparisons.items()
            ]
            if not np.isfinite(nonthermal_change):
                diagnostics.append(
                    {
                        "quantity": "nonthermal_rates",
                        "change": nonthermal_change,
                        "shell": None,
                        "state_index": None,
                        "current": np.nan,
                        "candidate": np.nan,
                        "absolute_change": np.nan,
                        "floor": None,
                        "unit": "comparison unavailable",
                    }
                )
            self.last_iteration_diagnostics = pd.DataFrame(
                diagnostics
            ).set_index("quantity")
            blocker = self.last_iteration_diagnostics["change"].idxmax()
            worst_entry = self.last_iteration_diagnostics.loc[blocker]
            blocker_report = (
                f"blocker={blocker}, shell={worst_entry['shell']}, "
                f"state_index={worst_entry['state_index']}, "
                f"current={worst_entry['current']:.6e}, "
                f"candidate={worst_entry['candidate']:.6e}, "
                f"absolute_change={worst_entry['absolute_change']:.6e} "
                f"{worst_entry['unit']}"
            )
            logger.info(
                "Nonthermal helium iteration %d/%d: maximum change %.3e; "
                "electron_density=%.3e, electron_temperature=%.3e, "
                "He II fraction delta=%.3e, level_population=%.3e "
                "(mixed, rtol=%.3e, atol=%.3e cm^-3), "
                "nonthermal_rates=%.3e; %s",
                iteration,
                self.settings.maximum_iterations,
                maximum_change,
                density_change,
                temperature_change,
                he_ii_fraction_change,
                level_change,
                self.settings.relative_tolerance,
                self.settings.absolute_population_tolerance,
                nonthermal_change,
                blocker_report,
            )
            active_factors = ion_factors.to_numpy()[
                plasma.number_density.loc[:, shell_columns].to_numpy() > 0.0
            ]
            if active_factors.size:
                factor_report = "unavailable"
                state_index = worst_entry["state_index"]
                if (
                    isinstance(state_index, tuple)
                    or blocker == "he_ii_fraction"
                ):
                    atomic_number = (
                        state_index[0] if isinstance(state_index, tuple) else 2
                    )
                    group = (atomic_number, worst_entry["shell"])
                    if group in ion_controller.diagnostics.index:
                        entry = ion_controller.diagnostics.loc[group]
                        factor_report = (
                            f"Z={group[0]}, shell={group[1]}, "
                            f"factor={entry['factor']:.3f}, reason={entry['reason']}"
                        )
                reasons = ion_controller.diagnostics["reason"]
                logger.info(
                    "Ion relaxation iteration %d: factor min/median/max="
                    "%.3f/%.3f/%.3f; full steps=%d/%d, promotions=%d, "
                    "backoffs=%d; blocker factor: %s",
                    iteration,
                    float(np.min(active_factors)),
                    float(np.median(active_factors)),
                    float(np.max(active_factors)),
                    int(np.count_nonzero(active_factors == 1.0)),
                    active_factors.size,
                    int((reasons == "promote").sum()),
                    int((reasons == "backoff").sum()),
                    factor_report,
                )
            if zero_target_count:
                largest_target = self.last_zero_target_levels.iloc[0]
                logger.info(
                    "Nonthermal helium level diagnostic: %d failing level/shell "
                    "entries; zero/below-atol targets=%d (atol=%.3e cm^-3); "
                    "largest absolute zero-target discrepancy: level=%s, "
                    "shell=%s, current=%.6e, candidate=%.6e, "
                    "absolute_change=%.6e cm^-3, allowed_change=%.6e cm^-3, "
                    "scaled_change=%.3e; "
                    "predicted relaxed population=%.6e cm^-3, "
                    "relaxed/current=%.6f (undamped level fractions, "
                    "damped ion totals).",
                    np.count_nonzero(failing_levels),
                    zero_target_count,
                    self.settings.absolute_population_tolerance,
                    largest_target["state_index"],
                    largest_target["shell"],
                    largest_target["current"],
                    largest_target["candidate"],
                    largest_target["absolute_change"],
                    largest_target["allowed_change"],
                    largest_target["scaled_change"],
                    largest_target["relaxed_population"],
                    largest_target["relaxed_to_current_ratio"],
                )
            terminal_nonthermal = nonthermal_state
            terminal_solution = solution
            if (
                other_maximum_change < self.settings.relative_tolerance
                and not np.any(failing_levels)
            ):
                current = candidate
                break
            current = relaxed_state
            previous_fixed_rates = nonthermal_state.fixed_rates
        else:
            worst = convergence_rows[-1]["maximum_change"]
            raise RuntimeError(
                "Nonthermal helium coupling did not converge after "
                f"{self.settings.maximum_iterations} iterations; "
                f"maximum change={worst:.3e}; {blocker_report}."
            )

        if terminal_nonthermal is None or terminal_solution is None:
            raise RuntimeError("No nonthermal helium iteration was performed.")
        plasma.update(
            electron_densities=pd.Series(
                current.electron_density, index=shell_columns
            ),
            ion_number_density=current.ion_number_density,
            level_number_density=current.level_number_density,
            link_t_rad_t_electron=(
                current.electron_temperature / complete_radiation_temperature
            ),
        )
        history = pd.DataFrame(convergence_rows).set_index("iteration")
        self.last_result = NonthermalCouplingResult(
            current,
            terminal_nonthermal,
            terminal_solution,
            history,
        )
        return self.last_result


class NonthermalHeliumOpticalWorkflow(StandardTARDISWorkflow):
    """Continuum optical workflow specialized for helium coupling."""

    plasma_coupling_solver: NonthermalHeliumPlasmaSolver

    def __init__(
        self,
        configuration: object,
        enable_virtual_packet_logging: bool = False,
        log_level: str | None = None,
        specific_log_level: str | None = None,
        show_progress_bars: bool = True,
        show_convergence_plots: bool = False,
        convergence_plots_kwargs: Mapping[str, object] | None = None,
        csvy: bool = False,
        continuum_interaction_species: Collection[str] = ("He I", "He II"),
        case_b_species: Collection[tuple[int, int]] = (),
    ) -> None:
        """Initialize continuum transport and helium equilibrium.

        The standard workflow initialization is reused for its simulation and
        convergence state. Continuum packets use direct macro-atom traversal.

        Parameters
        ----------
        configuration : object
            TARDIS configuration.
        enable_virtual_packet_logging : bool, optional
            Enable virtual-packet logging.
        log_level : str, optional
            Global workflow log level.
        specific_log_level : str, optional
            Component-specific workflow log level.
        show_progress_bars : bool, optional
            Display progress bars.
        show_convergence_plots : bool, optional
            Display convergence plots.
        convergence_plots_kwargs : mapping, optional
            Convergence-plot configuration.
        csvy : bool, optional
            Read the ejecta model from CSVY.
        continuum_interaction_species : collection of str, optional
            Species whose photoionization data enter continuum transport.
        case_b_species : collection of tuple[int, int], optional
            Species for which direct ground-state recombination is suppressed.
            The default is Case A for helium.
        """
        super().__init__(
            configuration,
            enable_virtual_packet_logging=enable_virtual_packet_logging,
            log_level=log_level,
            specific_log_level=specific_log_level,
            show_progress_bars=show_progress_bars,
            show_convergence_plots=show_convergence_plots,
            convergence_plots_kwargs=convergence_plots_kwargs,
            csvy=csvy,
        )
        self.configuration = configuration
        self.case_b_species = frozenset(case_b_species)
        montecarlo_globals.CONTINUUM_PROCESSES_ENABLED = True

        # Atomic data preparation is one-shot; continuum selection needs a
        # fresh owner after standard workflow initialization.
        self.atom_data = parse_atom_data(configuration)

        radiation_field = DilutePlanckianRadiationField(
            self.simulation_state.t_radiative,
            self.simulation_state.dilution_factor,
        )
        self.simulation_state.radiation_field_state = radiation_field
        elemental_number_density = (
            self.simulation_state.calculate_elemental_number_density(
                self.atom_data.atom_data.mass
            )
        )
        factory = PlasmaSolverFactory(self.atom_data, configuration)
        factory.helium_treatment = "none"
        factory.line_interaction_type = "macroatom"
        factory.continuum_interaction_species = list(
            continuum_interaction_species
        )
        factory.prepare_factory(
            self.simulation_state.abundance.index,
            "tardis.plasma.properties.property_collections",
            configuration,
            allow_continuum=True,
        )
        self.plasma_solver_factory = factory
        self.plasma_solver = factory.assemble(
            elemental_number_density,
            radiation_field,
            self.simulation_state.time_explosion,
            self.simulation_state._electron_densities,
        )
        self.plasma_solver.freeze(
            "electron_densities",
            "ion_number_density",
            "level_number_density",
        )

        self.opacity_solver = OpacitySolver(
            "macroatom", configuration.plasma.disable_line_scattering
        )
        self.macro_atom_solver = ContinuumMacroAtomSolver(
            self.atom_data.levels,
            self.atom_data.lines,
            self.atom_data.photoionization_data,
            self.atom_data.ionization_data,
            selected_continuum_transitions=np.asarray(
                factory.continuum_interaction_species_multi_index.tolist()
            ),
            line_interaction_type="macroatom",
            nthreads=configuration.plasma.continuum_interaction.nthreads,
            markov_chain_mode="explicit",
        )
        self.transport_state = None
        self.transport_solver = MCTransportSolverIIP.from_config(
            configuration,
            packet_source=self.simulation_state.packet_source,
            enable_virtual_packet_logging=enable_virtual_packet_logging,
        )
        self.transport_solver.line_interaction_type = "macroatom"
        self._continuum_estimators = None
        initial_state = self._current_plasma_state()
        maximum_electron_density = (
            self.plasma_solver.number_density.multiply(
                self.plasma_solver.number_density.index.to_numpy(dtype=float),
                axis=0,
            )
            .sum()
            .to_numpy(dtype=float)
        )
        bootstrap_evaluator = self._build_helium_evaluator(
            initial_state, None, (2, 0), maximum_electron_density
        )
        bootstrap_coefficients = (
            bootstrap_evaluator.calculate_continuum_coefficients(
                initial_state.electron_temperature
            )
        )
        self._build_continuum_states(
            bootstrap_coefficients[4], bootstrap_coefficients[1]
        )

    def _current_plasma_state(self) -> CoupledPlasmaState:
        """Return an immutable copy of the current continuum plasma state."""
        return CoupledPlasmaState(
            np.asarray(self.plasma_solver.electron_densities, dtype=np.float64),
            np.asarray(self.plasma_solver.t_electrons, dtype=np.float64),
            self.plasma_solver.ion_number_density.copy(deep=True),
            self.plasma_solver.level_number_density.copy(deep=True),
        )

    def _build_continuum_states(
        self,
        continuum_coefficients: ContinuumCoefficientState,
        level_to_continuum_saha_factor: pd.DataFrame,
    ) -> None:
        """Build continuum packet inputs from the accepted plasma state."""
        plasma = self.plasma_solver
        stimulated_cooling = None
        if self._continuum_estimators is not None:
            stimulated_cooling = pd.DataFrame(
                self._continuum_estimators.stim_recomb_cooling_estimator,
                index=plasma.level2continuum_idx.index,
                columns=plasma.number_density.columns,
            )
        self.continuum_macro_atom_state = (
            ContinuumMacroAtomState.from_equilibrium(
                plasma.atomic_data,
                plasma.lines,
                plasma.photo_ion_cross_sections,
                continuum_coefficients,
                level_to_continuum_saha_factor,
                plasma.ion_number_density,
                plasma.level_number_density,
                plasma.electron_densities,
                plasma.t_electrons,
                stimulated_cooling,
            )
        )
        self.continuum_opacity_state = ContinuumOpacityState.from_equilibrium(
            plasma.nu_i,
            plasma.level2continuum_idx,
            plasma.photo_ion_cross_sections,
            plasma.photo_ion_idx,
            plasma.ion_number_density,
            plasma.level_number_density,
            plasma.electron_densities,
            plasma.t_electrons,
            level_to_continuum_saha_factor,
            self.continuum_macro_atom_state.radiative_recombination_rate,
        )

    def _capture_continuum_estimators(self) -> None:
        """Normalize the fixed packet estimators for the next plasma solve."""
        transport = self.transport_state
        raw = transport.estimators_continuum
        volume = transport.geometry_state_numba.volume
        time = transport.time_of_simulation
        normalization = 1.0 / (time * volume * const.h.cgs.value).value
        target_intensity = (
            self.simulation_state.dilution_factor
            * self.simulation_state.t_radiative.value**4
            * const.sigma_sb.cgs.value
            / np.pi
        )
        estimated_intensity = transport.estimators_bulk.mean_intensity_total / (
            4.0 * np.pi * time.value * volume
        )
        damping = np.divide(
            target_intensity,
            estimated_intensity,
            out=np.ones_like(target_intensity),
            where=estimated_intensity > 0.0,
        )
        estimator = init_estimators_continuum(
            raw.photo_ion_estimator.shape,
            len(self.plasma_solver.number_density.columns),
        )
        estimator.photo_ion_estimator = (
            raw.photo_ion_estimator.copy()
            * normalization
            * damping
            * const.h.cgs.value
        )
        estimator.photo_ion_estimator_statistics = (
            raw.photo_ion_estimator_statistics.copy()
        )
        estimator.stim_recomb_estimator = (
            raw.stim_recomb_estimator.copy()
            * normalization
            * damping
            * const.h.cgs.value
        )
        estimator.bf_heating_estimator = (
            raw.bf_heating_estimator.copy()
            * normalization
            * const.h.cgs.value
            * damping
        )
        estimator.stim_recomb_cooling_estimator = (
            raw.stim_recomb_cooling_estimator.copy()
            * normalization
            * const.h.cgs.value
            * damping
        )
        ion_population = self.plasma_solver.ion_number_density
        charge_squared = np.square(
            ion_population.index.get_level_values("ion_number").to_numpy()
        )
        free_free_denominator = (
            self.plasma_solver.electron_densities.to_numpy()
            * ion_population.multiply(charge_squared, axis=0).sum().to_numpy()
        )
        free_free_factor = np.divide(
            np.sqrt(self.plasma_solver.t_electrons),
            free_free_denominator,
            out=np.zeros_like(free_free_denominator),
            where=free_free_denominator > 0.0,
        )
        estimator.ff_heating_estimator = (
            raw.ff_heating_estimator.copy()
            * free_free_factor
            * normalization
            * const.h.cgs.value
            * damping
        )
        self._continuum_estimators = estimator

        # Damp packet intensities before the estimator solver substitutes
        # Planck/epsilon fallbacks already evaluated in the target field.
        j_blues = self.transport_solver.radfield_prop_solver.estimate_jblues(
            transport.estimators_line.mean_intensity_blueward * damping,
            self.simulation_state.radiation_field_state,
            transport.time_explosion,
            time,
            volume,
            transport.opacity_state_numba.line_list_nu,
            detailed_optical_window=True,
        )
        self.plasma_solver.update(
            j_blues=pd.DataFrame(
                j_blues,
                index=self.plasma_solver.lines.index,
                columns=self.plasma_solver.number_density.columns,
            )
        )

    def _update_sobolev_state(self) -> None:
        """Refresh line optical depths after accepting a plasma state."""
        self._tau_sobolev = calculate_sobolev_line_opacity(
            self.atom_data.lines,
            self.plasma_solver.level_number_density,
            self.plasma_solver.time_explosion,
            self.plasma_solver.stimulated_emission_factor,
        )
        self._beta_sobolev = calculate_beta_sobolev(self._tau_sobolev)

    def _build_helium_evaluator(
        self,
        plasma_state: CoupledPlasmaState,
        fixed_rates: FixedRateState | None,
        bound_bound_species_filter: tuple[int, int] | None,
        maximum_electron_density: npt.NDArray[np.float64],
    ) -> PlasmaEquilibriumEvaluator:
        """Build one helium evaluator from the current adjacent-ion budget.

        Parameters
        ----------
        plasma_state : tardis.workflows.nonthermal_helium_workflow.CoupledPlasmaState
            Current shell-aligned ion and absolute level populations. He I
            plus He II sets the reduced population budget; He III is excluded.
        fixed_rates : tardis.plasma.equilibrium.inputs.FixedRateState, optional
            Frozen nonthermal excitation and ionization rates.
        bound_bound_species_filter : tuple[int, int], optional
            Atomic and ion number restricting bound-bound data, or all species
            when None.
        maximum_electron_density : numpy.ndarray
            Fully ionized electron-density bounds [cm^-3].

        Returns
        -------
        tardis.plasma.equilibrium.evaluator.PlasmaEquilibriumEvaluator
            Evaluator with fixed analytic or packet rates and current geometry.
        """
        plasma = self.plasma_solver
        shell_columns = plasma_state.ion_number_density.columns
        shell_positions = plasma.number_density.columns.get_indexer(
            shell_columns
        )
        if np.any(shell_positions < 0):
            raise ValueError(
                "Coupled shell labels are absent from the optical plasma."
            )
        equilibrium_species = (2, 0)
        level_index = plasma.level_number_density.index
        equilibrium_levels = plasma.atomic_data.levels.loc[level_index]
        photoionization_data = plasma.photo_ion_cross_sections.loc[
            plasma.photo_ion_cross_sections.index.isin(
                plasma.level2continuum_idx.index
            )
        ]
        level2continuum_edge_idx = plasma.level2continuum_idx
        species_level_positions = np.flatnonzero(
            (level_index.get_level_values("atomic_number") == 2)
            & (level_index.get_level_values("ion_number") == 0)
        )
        adjacent_ion_number_density = plasma_state.ion_number_density.loc[
            [
                equilibrium_species,
                (equilibrium_species[0], equilibrium_species[1] + 1),
            ]
        ].sum(axis=0)
        shell_number_densities = tuple(
            ShellNumberDensity(
                float(adjacent_ion_number_density.loc[shell]),
                plasma_state.level_number_density[shell].to_numpy(
                    dtype=np.float64
                ),
                species_level_positions,
            )
            for shell in shell_columns
        )

        all_line_index = plasma.lines.index
        all_line_species = all_line_index.droplevel(
            ["level_number_lower", "level_number_upper"]
        )
        if bound_bound_species_filter is None:
            line_positions = np.arange(len(all_line_index), dtype=np.int64)
        else:
            line_positions = np.flatnonzero(
                np.asarray(
                    all_line_species.isin([bound_bound_species_filter]),
                    dtype=bool,
                )
            )
        if len(line_positions) == 0:
            raise ValueError(
                "No bound-bound lines match "
                f"bound_bound_species_filter={bound_bound_species_filter}."
            )
        lines = plasma.lines.iloc[line_positions]
        line_index = lines.index
        lower_line_positions = plasma.lines_lower_level_index[line_positions]
        upper_line_positions = plasma.lines_upper_level_index[line_positions]
        sobolev_input = SobolevInputs(
            lower_line_positions,
            upper_line_positions,
            plasma.g.iloc[lower_line_positions].to_numpy(),
            plasma.g.iloc[upper_line_positions].to_numpy(),
            plasma.metastability.iloc[upper_line_positions].to_numpy(),
            np.asarray(
                line_index.droplevel(
                    ["level_number_lower", "level_number_upper"]
                ).isin([equilibrium_species]),
                dtype=bool,
            ),
            (
                lines.wavelength_cm.to_numpy()
                * lines.f_lu.to_numpy()
                * SOBOLEV_COEFFICIENT
                * plasma.time_explosion.to_value(u.s)
            ),
            np.arange(len(lines), dtype=np.int64),
            line_index,
        )
        collision_strengths = plasma.atomic_data.yg_data
        if bound_bound_species_filter is not None:
            collision_strengths = collision_strengths.loc[
                (
                    bound_bound_species_filter[0],
                    bound_bound_species_filter[1],
                    slice(None),
                    slice(None),
                ),
                :,
            ]
        rate_matrix_solver = RateMatrix(
            RadiativeRatesSolver(lines),
            ThermalCollisionalRateSolver(
                equilibrium_levels,
                lines,
                plasma.atomic_data.collision_data_temperatures,
                collision_strengths,
                collision_strengths_type="cmfgen",
            ),
            equilibrium_levels,
        )
        selected_number_density = plasma.number_density.loc[:, shell_columns]
        electron_density = pd.Series(
            plasma_state.electron_density, index=shell_columns
        )
        lte_ion_population, lte_level_population = calculate_lte_populations(
            plasma.thermal_phi_lte.loc[:, shell_columns],
            plasma.thermal_lte_partition_function.loc[:, shell_columns],
            selected_number_density,
            electron_density,
            plasma.thermal_lte_level_boltzmann_factor.loc[:, shell_columns],
            equilibrium_levels,
        )

        estimators = self._continuum_estimators
        if estimators is None:
            ion_rate_matrix = AnalyticIonRateMatrix(
                AnalyticPhotoionizationRateSolver(photoionization_data),
                CollisionalIonizationRateSolver(photoionization_data),
                fixed_rates=fixed_rates,
            )
        else:
            ion_rate_matrix = EstimatedIonRateMatrix(
                EstimatedPhotoionizationRateSolver(
                    photoionization_data,
                    level2continuum_edge_idx,
                    estimators,
                    1.0 * u.s,
                    1.0 * u.cm**3,
                ),
                CollisionalIonizationRateSolver(photoionization_data),
                plasma.phi.loc[:, shell_columns],
                fixed_rates=fixed_rates,
            )
        ion_population_solver = IonPopulationSolver(ion_rate_matrix)

        collisional_index = rate_matrix_solver.electron_rate_solver.all_collisional_strengths_index
        lower_level_index = pd.MultiIndex.from_arrays(
            [
                collisional_index.get_level_values("atomic_number"),
                collisional_index.get_level_values("ion_number"),
                collisional_index.get_level_values("level_number_lower"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        upper_level_index = pd.MultiIndex.from_arrays(
            [
                collisional_index.get_level_values("atomic_number"),
                collisional_index.get_level_values("ion_number"),
                collisional_index.get_level_values("level_number_upper"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        collisional_line_frequencies = (
            equilibrium_levels.energy.loc[upper_level_index].to_numpy()
            - equilibrium_levels.energy.loc[lower_level_index].to_numpy()
        ) / const.h.cgs.value
        radiation_field = DilutePlanckianRadiationField(
            plasma.dilute_planckian_radiation_field.temperature[
                shell_positions
            ],
            plasma.dilute_planckian_radiation_field.dilution_factor[
                shell_positions
            ],
        )
        thermal_balance_solver = ThermalBalanceSolver(
            BoundFreeThermalRates(photoionization_data),
            FreeFreeThermalRates(),
            CollisionalIonizationThermalRates(photoionization_data),
            CollisionalBoundThermalRates(
                pd.DataFrame({"nu": collisional_line_frequencies})
            ),
        )
        thermal_balance_arguments = {
            "free_free_heating_estimator": (
                pd.Series(0.0, index=shell_columns)
                if estimators is None
                else pd.Series(
                    estimators.ff_heating_estimator, index=shell_columns
                )
            ),
        }
        if estimators is None:
            thermal_balance_arguments["radiation_field"] = radiation_field
        else:
            thermal_balance_arguments["bound_free_heating_estimator"] = (
                pd.DataFrame(
                    estimators.bf_heating_estimator,
                    index=level2continuum_edge_idx.index,
                    columns=shell_columns,
                )
            )
            thermal_balance_arguments["stimulated_recombination_estimator"] = (
                pd.DataFrame(
                    estimators.stim_recomb_cooling_estimator,
                    index=level2continuum_edge_idx.index,
                    columns=shell_columns,
                )
            )
        return PlasmaEquilibriumEvaluator(
            photoionization_data,
            level2continuum_edge_idx,
            estimators,
            None if estimators is None else 1.0 * u.s,
            None if estimators is None else 1.0 * u.cm**3,
            equilibrium_levels,
            plasma.ionization_data,
            rate_matrix_solver,
            pd.DataFrame(
                np.asarray(plasma.j_blues)[line_positions][:, shell_positions],
                index=line_index,
                columns=shell_columns,
            ),
            shell_number_densities,
            tuple(sobolev_input for _ in shell_columns),
            level_index,
            equilibrium_species,
            selected_number_density,
            maximum_electron_density,
            ion_population_solver=ion_population_solver,
            ion_population_arguments={
                "radiation_field": (
                    radiation_field if estimators is None else None
                ),
                "elemental_number_density": selected_number_density,
                "lte_level_population": lte_level_population,
                "lte_ion_population": lte_ion_population,
                "estimated_ion_population": plasma_state.ion_number_density,
                "partition_function": plasma.thermal_lte_partition_function.loc[
                    :, shell_columns
                ],
                "boltzmann_factor": plasma.thermal_lte_level_boltzmann_factor.loc[
                    :, shell_columns
                ],
            },
            thermal_balance_solver=thermal_balance_solver,
            thermal_balance_arguments=thermal_balance_arguments,
            radiation_field=(radiation_field if estimators is None else None),
            fixed_rates=fixed_rates,
            case_b_species=self.case_b_species,
            time_explosion=plasma.time_explosion,
            dilute_radiation=DiluteRadiationInputs(
                np.asarray(plasma.t_rad, dtype=np.float64)[
                    shell_positions
                ].copy(),
                np.asarray(plasma.w, dtype=np.float64)[shell_positions].copy(),
                plasma.zeta_data,
                plasma.level_boltzmann_factor.loc[:, shell_columns].copy(
                    deep=True
                ),
                plasma.partition_function.loc[:, shell_columns].copy(deep=True),
                self.configuration.plasma.get("delta_treatment", None),
            ),
            background_level_population_fraction=(
                plasma.level_boltzmann_factor.loc[level_index, shell_columns]
                / reindex_ion_population_to_level_population(
                    plasma.partition_function.loc[:, shell_columns],
                    plasma.level_boltzmann_factor.loc[
                        level_index, shell_columns
                    ],
                    next_higher=False,
                )
            ),
        )

    def build_helium_evaluator(
        self,
        plasma_state: CoupledPlasmaState,
        fixed_rates: FixedRateState,
        bound_bound_species_filter: tuple[int, int] | None,
    ) -> PlasmaEquilibriumEvaluator:
        """Build a fixed estimator evaluator for one Picard iteration."""
        maximum_electron_density = (
            self.plasma_solver.number_density.multiply(
                self.plasma_solver.number_density.index.to_numpy(dtype=float),
                axis=0,
            )
            .sum()
            .to_numpy(dtype=float)
        )
        return self._build_helium_evaluator(
            plasma_state,
            fixed_rates,
            bound_bound_species_filter,
            maximum_electron_density,
        )

    def solve_plasma(self, estimated_radfield_properties: object) -> None:
        """Freeze continuum estimators, then converge helium."""
        super().solve_plasma(estimated_radfield_properties)
        self._capture_continuum_estimators()
        result = self.plasma_coupling_solver.solve(self)
        evaluation = result.equilibrium_solution.evaluation
        self._build_continuum_states(
            evaluation.continuum_coefficients,
            evaluation.level_to_continuum_saha_factor,
        )
        self._update_sobolev_state()

    def solve_opacity(self) -> dict[str, object]:
        """Solve line and continuum opacity and direct macro-atom states."""
        self._update_sobolev_state()
        tau_sobolev = self._tau_sobolev
        beta_sobolev = self._beta_sobolev
        if self.opacity_solver.disable_line_scattering:
            tau_sobolev = pd.DataFrame(
                np.zeros_like(tau_sobolev),
                index=tau_sobolev.index,
                columns=tau_sobolev.columns,
            )
            beta_sobolev = pd.DataFrame(
                np.ones_like(beta_sobolev),
                index=beta_sobolev.index,
                columns=beta_sobolev.columns,
            )
        opacity_state = self.opacity_solver.solve(
            self.plasma_solver,
            tau_sobolev=tau_sobolev,
            beta_sobolev=beta_sobolev,
        )
        macro_atom_state = self.macro_atom_solver.solve(
            self.plasma_solver.j_blues,
            opacity_state.beta_sobolev,
            self.plasma_solver.stimulated_emission_factor,
            self.continuum_macro_atom_state,
            self.plasma_solver.electron_densities,
        )
        return {
            "opacity_state": opacity_state,
            "macro_atom_state": macro_atom_state,
        }

    def solve_montecarlo(
        self,
        opacity_states: dict[str, object],
        no_of_real_packets: int,
        no_of_virtual_packets: int = 0,
    ) -> None:
        """Run continuum Monte Carlo transport and collect estimators."""
        self.transport_state = self.transport_solver.initialize_transport_state(
            self.simulation_state,
            opacity_states["opacity_state"],
            opacity_states["macro_atom_state"],
            self.continuum_opacity_state,
            no_of_real_packets,
            no_of_virtual_packets=no_of_virtual_packets,
            iteration=self.completed_iterations,
        )
        self.transport_solver.run(
            self.transport_state,
            show_progress_bars=self.show_progress_bars,
        )
        output_energy = self.transport_state.packet_collection.output_energies
        if np.sum(output_energy < 0) == len(output_energy):
            logger.critical("No r-packet escaped through the outer boundary.")


class NonthermalHeliumWorkflow:
    """Own high-energy and coupled optical nonthermal-helium calculations."""

    def __init__(
        self,
        optical_workflow: NonthermalHeliumOpticalWorkflow,
        high_energy_workflow: TARDISHEWorkflow,
        plasma_coupling_solver: NonthermalHeliumPlasmaSolver,
        high_energy_settings: HighEnergySettings,
    ) -> None:
        """Compose the optical, high-energy, and plasma-coupling owners."""
        self.optical_workflow = optical_workflow
        self.high_energy_workflow = high_energy_workflow
        self.plasma_coupling_solver = plasma_coupling_solver
        self.high_energy_settings = high_energy_settings
        self.high_energy_result: TARDISHEWorkflowResult | None = None
        self.deposition_state: DepositionState | None = None

    def run(self) -> None:
        """Run high-energy transport once, then the coupled optical workflow."""
        optical_state = self.optical_workflow.simulation_state
        high_energy_state = self.high_energy_workflow.simulation_state
        for boundary_name in ("v_inner", "v_outer"):
            optical_boundary = getattr(optical_state, boundary_name).to_value(
                u.cm / u.s
            )
            high_energy_boundary = getattr(
                high_energy_state, boundary_name
            ).to_value(u.cm / u.s)
            if (
                optical_boundary.shape != high_energy_boundary.shape
                or not np.allclose(
                    optical_boundary,
                    high_energy_boundary,
                    rtol=1.0e-12,
                    atol=0.0,
                )
            ):
                raise ValueError(
                    "High-energy and optical shell velocity boundaries do "
                    f"not match for {boundary_name}."
                )
        self.high_energy_result = self.high_energy_workflow.run(
            **self.high_energy_settings.as_run_arguments()
        )
        time_bin_edges, effective_times = self.high_energy_workflow.get_times(
            self.high_energy_settings.time_start,
            self.high_energy_settings.time_end,
            self.high_energy_settings.time_space,
            self.high_energy_settings.time_steps,
        )
        time_bin_edges = time_bin_edges * u.d
        effective_times = (effective_times * u.d).to_value(u.s)
        velocity_volume = calculate_ejecta_velocity_volume(
            self.high_energy_workflow.simulation_state
        )
        shell_volume = (
            velocity_volume[:, np.newaxis]
            * effective_times[np.newaxis, :] ** 3
            * u.cm**3
        )
        shell_index = self.optical_workflow.plasma_solver.number_density.columns
        self.deposition_state = prepare_deposition_state(
            self.high_energy_result,
            time_bin_edges,
            optical_state.time_explosion,
            shell_volume,
            shell_index,
        )
        self.plasma_coupling_solver.deposition_state = self.deposition_state
        self.optical_workflow.plasma_coupling_solver = (
            self.plasma_coupling_solver
        )
        self.optical_workflow.run()
