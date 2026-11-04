from pathlib import Path
from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u

from tardis import constants as const
from tardis.io.configuration.config_reader import Configuration
from tardis.plasma.equilibrium.inputs import FixedRateState
from tardis.plasma.equilibrium.solver import PlasmaEquilibriumSolution
from tardis.plasma.nonthermal.pynonthermal_solver import (
    NonthermalPlasmaState,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.transport.montecarlo.estimators import init_estimators_continuum
from tardis.transport.montecarlo.estimators.mc_rad_field_solver import (
    MCRadiationFieldPropertiesSolver,
)
from tardis.workflows.high_energy.deposition import DepositionState
from tardis.workflows.nonthermal_helium_workflow import (
    CoupledPlasmaState,
    HighEnergySettings,
    NonthermalCouplingSettings,
    NonthermalHeliumOpticalWorkflow,
    NonthermalHeliumPlasmaSolver,
    NonthermalHeliumWorkflow,
    _IonRelaxationController,
)
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow


class FakeNonthermalSolver:
    """Return rates that track the supplied electron density."""

    def __init__(self) -> None:
        self.electron_densities: list[float] = []
        self.electron_temperatures: list[np.ndarray] = []

    def solve(self, plasma_input: object) -> NonthermalPlasmaState:
        electron_densities = plasma_input.electron_number_density.to_value(
            u.cm**-3
        )
        electron_density = float(electron_densities[0])
        self.electron_densities.append(electron_density)
        self.electron_temperatures.append(
            plasma_input.electron_temperature.to_value(u.K).copy()
        )
        shell_columns = plasma_input.ion_number_density.columns
        ion_index = pd.MultiIndex.from_tuples(
            [(2, 0, 0, 1)],
            names=[
                "atomic_number",
                "ion_number",
                "ion_number_source",
                "ion_number_destination",
            ],
        )
        return NonthermalPlasmaState(
            FixedRateState(
                ionization_rates=pd.DataFrame(
                    [electron_densities / 1.0e9],
                    index=ion_index,
                    columns=shell_columns,
                ),
                heating_rate_density=pd.Series(
                    electron_densities / 1.0e9, index=shell_columns
                ),
            ),
            pd.DataFrame(columns=shell_columns),
        )


class FixedEquilibriumSolver:
    """Return the same known physical equilibrium on every call."""

    def __init__(self) -> None:
        self.level_seeds: list[pd.DataFrame] = []
        self.fixed_temperature_ratio: float | None = None

    def solve(
        self, *args: object, **kwargs: object
    ) -> PlasmaEquilibriumSolution:
        del kwargs
        level_seed = args[3]
        self.level_seeds.append(level_seed)
        shell_columns = level_seed.columns
        shell_count = len(shell_columns)
        ion_index = pd.MultiIndex.from_tuples(
            [(2, 0), (2, 1)], names=["atomic_number", "ion_number"]
        )
        level_index = pd.MultiIndex.from_tuples(
            [(2, 0, 0), (2, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        )
        evaluation = SimpleNamespace(
            ion_population=pd.DataFrame(
                np.tile([[6.0e8], [4.0e8]], shell_count),
                index=ion_index,
                columns=shell_columns,
            ),
            absolute_level_population=pd.DataFrame(
                np.tile([[4.0e8], [2.0e8]], shell_count),
                index=level_index,
                columns=shell_columns,
            ),
        )
        return PlasmaEquilibriumSolution(
            np.full(shell_count, 4.0e8),
            np.full(shell_count, 8000.0),
            evaluation,
        )


class FakePlasma:
    """Minimal mutable optical plasma used by the coupling test."""

    def __init__(self, shell_count: int = 1) -> None:
        shell_columns = pd.RangeIndex(shell_count)
        self.electron_densities = pd.Series(
            np.arange(shell_count) * 1.0e8 + 2.0e8, index=shell_columns
        )
        self.t_electrons = np.arange(shell_count) * 1000.0 + 6000.0
        ion_index = pd.MultiIndex.from_tuples(
            [(2, 0), (2, 1)], names=["atomic_number", "ion_number"]
        )
        self.ion_number_density = pd.DataFrame(
            np.tile([[8.0e8], [2.0e8]], shell_count),
            index=ion_index,
            columns=shell_columns,
        )
        level_index = pd.MultiIndex.from_tuples(
            [(2, 0, 0), (2, 0, 1)],
            names=["atomic_number", "ion_number", "level_number"],
        )
        self.level_number_density = pd.DataFrame(
            np.tile([[6.0e8], [2.0e8]], shell_count),
            index=level_index,
            columns=shell_columns,
        )
        self.number_density = pd.DataFrame(
            np.tile([[0.0], [1.0e9]], shell_count),
            index=[1, 2],
            columns=shell_columns,
        )

    def update(self, **kwargs: object) -> None:
        for name, value in kwargs.items():
            setattr(self, name, value)


@pytest.mark.parametrize(
    ("relaxation", "represented_candidate"),
    [(0.5, False), (0.75, False), (1.0, False), (0.25, True)],
)
@pytest.mark.parametrize("converged_input", [False, True])
@pytest.mark.parametrize("initial_parent", [0.0, np.nextafter(0.0, 1.0)])
def test_coupling_canonicalizes_underflowed_nickel_parent(
    relaxation: float,
    converged_input: bool,
    initial_parent: float,
    represented_candidate: bool,
) -> None:
    """An unrepresentable Ni IX distribution stays zero across coupling updates.

    Regime: all 26 Ni IX levels and the old parent are zero; the evaluator
    returns the smallest binary64 subnormal parent but zero levels, shell 36.
    The exact zero limiting state is the oracle, not a population tolerance.
    A second case returns a represented subnormal parent of 26 quanta with
    one quantum in every level: it must survive direct acceptance, but a
    quarter-step reconstruction underflows all its levels and must repair
    the relaxed parent. Both next-input and direct converged-commit paths
    must repair parent and
    levels without mutating evaluator output or changing representable nuclei
    or charge sums in a gas dominated by He and Ni VIII.
    """
    plasma = FakePlasma()
    for frame in (
        plasma.ion_number_density,
        plasma.level_number_density,
        plasma.number_density,
    ):
        frame.columns = pd.Index([36])
    plasma.electron_densities.index = pd.Index([36])
    plasma.number_density.loc[28] = 1.0e4
    plasma.ion_number_density.loc[(28, 7), :] = 1.0e4
    plasma.ion_number_density.loc[(28, 8), :] = initial_parent
    plasma.level_number_density.loc[(28, 7, 0), :] = 1.0e4
    for level_number in range(26):
        plasma.level_number_density.loc[(28, 8, level_number), :] = 0.0
    if converged_input:
        plasma.electron_densities[:] = 4.0e8
        plasma.t_electrons[:] = 8000.0
        plasma.ion_number_density.loc[(2, 0), :] = 6.0e8
        plasma.ion_number_density.loc[(2, 1), :] = 4.0e8
        plasma.level_number_density.loc[(2, 0, 0), :] = 4.0e8
    states: list[CoupledPlasmaState] = []
    evaluations: list[SimpleNamespace] = []

    class UnderflowEquilibriumSolver(FixedEquilibriumSolver):
        def solve(
            self, *args: object, **kwargs: object
        ) -> PlasmaEquilibriumSolution:
            solution = super().solve(*args, **kwargs)
            evaluation = solution.evaluation
            evaluation.ion_population.loc[(28, 7), :] = 1.0e4
            quantum = np.nextafter(0.0, 1.0)
            evaluation.ion_population.loc[(28, 8), :] = (
                26 * quantum if represented_candidate else quantum
            )
            evaluation.absolute_level_population.loc[(28, 7, 0), :] = 1.0e4
            for level_number in range(26):
                evaluation.absolute_level_population.loc[
                    (28, 8, level_number), :
                ] = quantum if represented_candidate else 0.0
            evaluations.append(evaluation)
            return solution

    def evaluator_factory(
        state: CoupledPlasmaState, *args: object
    ) -> SimpleNamespace:
        del args
        states.append(state)
        return SimpleNamespace()

    coupling = NonthermalHeliumPlasmaSolver(
        FakeNonthermalSolver(),
        evaluator_factory,
        equilibrium_solver=UnderflowEquilibriumSolver(),
        settings=NonthermalCouplingSettings(
            relaxation=relaxation,
            adaptive_ion_relaxation=False,
            maximum_iterations=60,
        ),
    )
    coupling.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        pd.Index([36]),
        np.ones(1) * u.erg / u.s / u.cm**3,
        np.zeros(1) * u.erg / u.s / u.cm**3,
        np.ones(1) * u.erg / u.s / u.cm**3,
    )
    result = coupling.solve(
        SimpleNamespace(
            plasma_solver=plasma,
            simulation_state=SimpleNamespace(t_radiative=np.array([1.0e4])),
        )
    )
    # Fixed-rate comparison requires a previous iteration even at equilibrium.
    assert len(states) == (2 if converged_input else len(evaluations))
    quantum = np.nextafter(0.0, 1.0)
    for state in [*states, coupling.last_candidate_state, result.plasma_state]:
        is_candidate = (
            state is coupling.last_candidate_state
            or state is result.plasma_state
        )
        expected_parent = (
            26 * quantum if represented_candidate and is_candidate else 0.0
        )
        npt.assert_equal(
            state.ion_number_density.loc[(28, 8), 36], expected_parent
        )
        npt.assert_array_equal(
            state.level_number_density.loc[(28, 8), 36],
            quantum if represented_candidate and is_candidate else 0.0,
        )
        nickel = state.ion_number_density.loc[28, 36]
        npt.assert_equal(nickel.sum(), 1.0e4)
        npt.assert_equal(nickel.multiply(nickel.index).sum(), 7.0e4)
    npt.assert_equal(
        plasma.ion_number_density.loc[(28, 8), 36],
        26 * quantum if represented_candidate else 0.0,
    )
    for evaluation in evaluations:
        npt.assert_equal(
            evaluation.ion_population.loc[(28, 8), 36],
            26 * quantum if represented_candidate else quantum,
        )


@pytest.mark.parametrize("stage", ["initial", "candidate", "relaxed"])
def test_nonthermal_rejects_nonfinite_levels_at_first_boundary(
    stage: str,
) -> None:
    """Report the first invalid level without another degradation solve.

    A NaN in He I level 1, shell 0, or an undefined He II level fraction
    introduced by relaxation must not be treated as zero or committed. This
    is orchestration/error-localization verification; no physical population
    solution is inferred from the synthetic evaluator.
    """
    provider = FakeNonthermalSolver()
    plasma = FakePlasma()
    if stage == "initial":
        plasma.level_number_density.loc[(2, 0, 1), 0] = np.nan
    if stage == "relaxed":
        plasma.level_number_density.loc[(2, 1, 0), :] = 0.0
    initial = plasma.level_number_density.copy(deep=True)

    class NonfiniteEquilibriumSolver(FixedEquilibriumSolver):
        def solve(
            self, *args: object, **kwargs: object
        ) -> PlasmaEquilibriumSolution:
            solution = super().solve(*args, **kwargs)
            if stage == "relaxed":
                solution.evaluation.absolute_level_population.loc[
                    (2, 1, 0), :
                ] = 0.0
            else:
                solution.evaluation.absolute_level_population.loc[
                    (2, 0, 1), 0
                ] = np.nan
            return solution

    def evaluator_factory(*args: object) -> SimpleNamespace:
        del args
        return SimpleNamespace()

    coupling = NonthermalHeliumPlasmaSolver(
        provider,
        evaluator_factory,
        equilibrium_solver=NonfiniteEquilibriumSolver(),
    )
    coupling.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        pd.Index([0]),
        np.ones(1) * u.erg / u.s / u.cm**3,
        np.zeros(1) * u.erg / u.s / u.cm**3,
        np.ones(1) * u.erg / u.s / u.cm**3,
    )
    with pytest.raises(
        ValueError,
        match=(
            rf"Nonfinite {stage} level population.*level="
            + (r"\(2, 1, 0\)" if stage == "relaxed" else r"\(2, 0, 1\)")
            + r".*shell=0"
        ),
    ):
        coupling.solve(
            SimpleNamespace(
                plasma_solver=plasma,
                simulation_state=SimpleNamespace(t_radiative=np.array([1.0e4])),
            )
        )
    assert len(provider.electron_densities) == (0 if stage == "initial" else 1)
    assert coupling.last_result is None
    pdt.assert_frame_equal(plasma.level_number_density, initial)


@pytest.mark.parametrize("adaptive", [True, False])
def test_adaptive_ion_relaxation_removes_constant_candidate_tail(
    adaptive: bool,
) -> None:
    """A stable candidate removes a trace-ion tail without changing ne damping.

    Regime: a synthetic Li population moves entirely from Li I to Li II,
    while the helium/electron root remains fixed. With alpha=0.5 the Li I
    tail needs more than twenty steps to reach atol=1e-20 cm^-3. Verification:
    adaptive factors [0.5, 0.5, 0.75, 0.75, 1] give exactly zero Li I in the
    sixth input, while ne follows 4e8 - 2e8 * 0.5**k. Comparisons use exact
    zero or rtol=1e-14 for roundoff. This verifies fixed-point orchestration,
    not atomic rates or model performance.
    """
    plasma = FakePlasma()
    plasma.number_density.loc[3] = 1.0
    plasma.ion_number_density.loc[(3, 0), :] = 1.0
    plasma.ion_number_density.loc[(3, 1), :] = 0.0
    plasma.level_number_density.loc[(3, 0, 0), :] = 1.0
    plasma.level_number_density.loc[(3, 1, 0), :] = 0.0
    states = []

    class ConstantCandidateSolver(FixedEquilibriumSolver):
        def solve(
            self, *args: object, **kwargs: object
        ) -> PlasmaEquilibriumSolution:
            solution = super().solve(*args, **kwargs)
            solution.evaluation.ion_population.loc[(3, 0), :] = 0.0
            solution.evaluation.ion_population.loc[(3, 1), :] = 1.0
            solution.evaluation.absolute_level_population.loc[(3, 0, 0), :] = (
                0.0
            )
            solution.evaluation.absolute_level_population.loc[(3, 1, 0), :] = (
                1.0
            )
            return solution

    def evaluator_factory(
        state: CoupledPlasmaState, *args: object
    ) -> SimpleNamespace:
        del args
        states.append(state)
        return SimpleNamespace()

    coupling = NonthermalHeliumPlasmaSolver(
        FakeNonthermalSolver(),
        evaluator_factory,
        equilibrium_solver=ConstantCandidateSolver(),
        settings=NonthermalCouplingSettings(adaptive_ion_relaxation=adaptive),
    )
    coupling.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        pd.Index([0]),
        np.ones(1) * u.erg / u.s / u.cm**3,
        np.zeros(1) * u.erg / u.s / u.cm**3,
        np.ones(1) * u.erg / u.s / u.cm**3,
    )
    optical = SimpleNamespace(
        plasma_solver=plasma,
        simulation_state=SimpleNamespace(t_radiative=np.array([1.0e4])),
    )
    if adaptive:
        result = coupling.solve(optical)
        assert len(result.convergence_history) == 10
        npt.assert_array_equal(states[5].ion_number_density.loc[(3, 0)], 0.0)
        npt.assert_array_equal(
            states[5].level_number_density.loc[(3, 0, 0)], 0.0
        )
    else:
        with pytest.raises(RuntimeError, match="did not converge after 20"):
            coupling.solve(optical)
        npt.assert_allclose(states[-1].ion_number_density.loc[(3, 0)], 0.5**19)
        assert coupling.last_result is None
    history = coupling.last_ion_relaxation_history.xs(
        (3, 0), level=["atomic_number", "shell"]
    )
    npt.assert_array_equal(
        history["factor"].iloc[:5],
        [0.5, 0.5, 0.75, 0.75, 1.0] if adaptive else [0.5] * 5,
    )
    for state in states:
        npt.assert_allclose(
            state.ion_number_density.groupby(level="atomic_number").sum(),
            plasma.number_density.loc[[2, 3]],
            rtol=1e-14,
        )
        assert (state.ion_number_density.to_numpy() >= 0.0).all()
        level_sums = state.level_number_density.groupby(
            level=["atomic_number", "ion_number"]
        ).sum()
        pdt.assert_frame_equal(
            level_sums,
            state.ion_number_density.reindex(level_sums.index),
            rtol=1e-14,
            atol=0.0,
        )
    npt.assert_allclose(
        [state.electron_density[0] for state in states[:6]],
        4.0e8 - 2.0e8 * 0.5 ** np.arange(6),
        rtol=1e-14,
    )
    if adaptive:
        coupling.solve(optical)
        first = coupling.last_ion_relaxation_history.loc[(1, 3, 0)]
        npt.assert_equal(first["factor"], 0.5)
        assert first["reason"] == "initial"


def test_ion_controller_promotes_and_backs_off_independently() -> None:
    """Verify exact controller factors from normalized candidate trajectories.

    Regime: two elements with abundances 1e9 and 1e-12 cm^-3, two reordered
    shells, and an absent element. Stable targets promote on comparisons 2
    and 4; material reversals and >2x growth back off to 0.25. Small reversals
    of 1e-4 in elemental fraction must not trigger backoff. Factors/counts are
    exact; abundance conservation allows rtol=1e-14 for floating-point sums.
    This numerical verification detects coupling controllers across labels,
    using raw densities, reacting to noise, and losing convex conservation.
    """
    densities = pd.DataFrame(
        [[1e-12, 2e-12], [1e9, 2e9], [0.0, 0.0]],
        index=[3, 2, 1],
        columns=[53, 25],
    )
    ion_index = pd.MultiIndex.from_tuples(
        [(3, 1), (2, 0), (1, 0), (2, 1), (3, 0), (1, 1)],
        names=["atomic_number", "ion_number"],
    )
    controller = _IonRelaxationController(
        NonthermalCouplingSettings(), densities
    )
    paths = {
        (2, 53): [0.4] * 5,
        (2, 25): [0.2, 0.4, 0.2, 0.4, 0.2],
        (3, 53): [0.1, 0.11, 0.14, 0.15, 0.16],
        (3, 25): [0.5, 0.5001, 0.5, 0.5001, 0.5],
    }
    current = pd.DataFrame(0.0, index=ion_index, columns=densities.columns)
    for atomic_number in [1, 2, 3]:
        current.loc[(atomic_number, 0)] = densities.loc[atomic_number]
    history = []
    for iteration in range(5):
        candidate = pd.DataFrame(
            0.0, index=ion_index, columns=densities.columns
        )
        for (atomic_number, shell), path in paths.items():
            total = densities.loc[atomic_number, shell]
            candidate.loc[(atomic_number, 0), shell] = total * path[iteration]
            candidate.loc[(atomic_number, 1), shell] = total * (
                1 - path[iteration]
            )
        original = candidate.copy(deep=True)
        factors = controller.update(candidate)
        history.append(controller.diagnostics.copy())
        stage_factors = pd.DataFrame(
            factors.reindex(
                ion_index.get_level_values("atomic_number")
            ).to_numpy(),
            index=ion_index,
            columns=densities.columns,
        )
        current = stage_factors * candidate + (1 - stage_factors) * current
        npt.assert_allclose(
            current.groupby(level="atomic_number").sum(),
            densities.sort_index(),
            rtol=1e-14,
            atol=0.0,
        )
        assert (current.to_numpy() >= 0).all()
        pdt.assert_frame_equal(candidate, original)
    npt.assert_array_equal(factors.loc[[2, 3]], [[1.0, 0.25], [0.25, 1.0]])
    assert history[2].loc[(2, 25), "reason"] == "backoff"
    assert history[2].loc[(3, 53), "reason"] == "backoff"
    assert history[2].loc[(3, 25), "reason"] == "promote"
    npt.assert_array_equal(factors.loc[1], 1.0)
    npt.assert_array_equal(
        controller.diagnostics.loc[1, "stable_comparisons"], 0
    )


@pytest.mark.parametrize(
    "settings",
    [
        {"minimum_ion_relaxation": 0.0},
        {"maximum_ion_relaxation": 1.1},
        {"ion_candidate_stability_tolerance": np.nan},
        {"ion_stable_comparisons": 0},
    ],
)
def test_adaptive_ion_settings_reject_invalid_bounds(
    settings: dict[str, float],
) -> None:
    """Prevent extrapolation factors and undefined stability criteria."""
    with pytest.raises(ValueError, match="Adaptive ion relaxation requires"):
        NonthermalCouplingSettings(**settings)


def test_nonthermal_rates_refresh_only_between_equilibrium_solves() -> None:
    """Refresh pynonthermal across Picard steps and commit convergence."""
    provider = FakeNonthermalSolver()
    equilibrium_solver = FixedEquilibriumSolver()
    factory_species_filters: list[tuple[int, int] | None] = []
    optical_workflow = SimpleNamespace(
        plasma_solver=FakePlasma(),
        simulation_state=SimpleNamespace(t_radiative=np.array([1.0e4])),
    )

    def evaluator_factory(
        state: object,
        fixed_rates: object,
        workflow: object,
        bound_bound_species_filter: tuple[int, int] | None,
    ) -> SimpleNamespace:
        del state, fixed_rates, workflow
        factory_species_filters.append(bound_bound_species_filter)
        return SimpleNamespace()

    coupling_solver = NonthermalHeliumPlasmaSolver(
        provider,
        evaluator_factory,
        equilibrium_solver=equilibrium_solver,
        bound_bound_species_filter=(2, 0),
        settings=NonthermalCouplingSettings(
            maximum_iterations=4,
            relative_tolerance=1.0e-12,
            relaxation=1.0,
        ),
    )
    coupling_solver.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        pd.Index([0]),
        np.array([1.0]) * u.erg / u.s / u.cm**3,
        np.array([0.0]) * u.erg / u.s / u.cm**3,
        np.array([1.0]) * u.erg / u.s / u.cm**3,
    )

    result = coupling_solver.solve(optical_workflow)

    npt.assert_array_equal(provider.electron_densities, [2.0e8, 4.0e8, 4.0e8])
    npt.assert_array_equal(result.plasma_state.electron_density, [4.0e8])
    npt.assert_array_equal(
        optical_workflow.plasma_solver.electron_densities, [4.0e8]
    )
    assert len(result.convergence_history) == 3
    # Nonzero targets must not be mistaken for zero-target relaxation stalls.
    npt.assert_array_equal(
        result.convergence_history["zero_target_level_count"], 0
    )
    assert coupling_solver.last_zero_target_levels.empty
    assert factory_species_filters == [(2, 0)] * 3
    for level_seed in equilibrium_solver.level_seeds:
        assert level_seed.index.equals(pd.Index([0, 1], name="level_number"))
        npt.assert_allclose(level_seed.sum(axis=0), 1.0)


def test_nonthermal_default_prescribes_temperature_before_rate_refresh() -> (
    None
):
    """Use 0.9 T_rad in pynonthermal, evaluator, relaxation, and plasma commit.

    Claim: prescribing temperature applies before the first degradation solve
    and remains exact throughout Picard relaxation. Regime: two shells whose
    initial temperatures and synthetic thermal root disagree with 0.9 T_rad.
    Verification: the recorded provider inputs and committed temperatures equal
    [9000, 10800] K; the independent linear charge root is 4e8 cm^-3.
    """
    provider = FakeNonthermalSolver()
    plasma = FakePlasma(shell_count=2)
    optical_workflow = SimpleNamespace(
        plasma_solver=plasma,
        simulation_state=SimpleNamespace(t_radiative=np.array([1.0e4, 1.2e4])),
    )
    evaluators = []

    class Evaluator:
        thermal_balance_solver = object()

        def evaluate(
            self,
            electron_density: np.ndarray,
            electron_temperature: np.ndarray,
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            return SimpleNamespace(
                electron_residual=pd.Series(
                    electron_density / 4.0e8 - 1.0, index=level_seed.columns
                ),
                fractional_heating=(
                    None
                    if self.thermal_balance_solver is None
                    else pd.Series(
                        electron_temperature / 8000.0 - 1.0,
                        index=level_seed.columns,
                    )
                ),
                normalized_population=level_seed,
                ion_population=pd.DataFrame(
                    [[6.0e8, 6.0e8], [4.0e8, 4.0e8]],
                    index=plasma.ion_number_density.index,
                    columns=level_seed.columns,
                ),
                absolute_level_population=pd.DataFrame(
                    [[4.0e8, 4.0e8], [2.0e8, 2.0e8]],
                    index=plasma.level_number_density.index,
                    columns=level_seed.columns,
                ),
            )

    def evaluator_factory(
        *args: object,
    ) -> Evaluator:
        del args
        evaluator = Evaluator()
        evaluators.append(evaluator)
        return evaluator

    coupling_solver = NonthermalHeliumPlasmaSolver(
        provider,
        evaluator_factory,
        settings=NonthermalCouplingSettings(
            maximum_iterations=40, relaxation=0.5
        ),
    )
    coupling_solver.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        plasma.number_density.columns,
        np.ones(2) * u.erg / u.s / u.cm**3,
        np.zeros(2) * u.erg / u.s / u.cm**3,
        np.ones(2) * u.erg / u.s / u.cm**3,
    )

    result = coupling_solver.solve(optical_workflow)

    for temperature in provider.electron_temperatures:
        npt.assert_array_equal(temperature, [9000.0, 10800.0])
    npt.assert_array_equal(
        result.plasma_state.electron_temperature, [9000.0, 10800.0]
    )
    npt.assert_array_equal(plasma.link_t_rad_t_electron, [0.9, 0.9])
    npt.assert_allclose(result.plasma_state.electron_density, 4.0e8, rtol=1e-7)
    assert all(
        evaluator.thermal_balance_solver is None for evaluator in evaluators
    )
    assert len(provider.electron_temperatures) > 1
    assert (
        (result.nonthermal_state.fixed_rates.ionization_rates > 0.0).all().all()
    )

    # The next optical iteration supplies a different radiation temperature.
    previous_call_count = len(provider.electron_temperatures)
    optical_workflow.simulation_state.t_radiative = np.array([8000.0, 11000.0])
    next_result = coupling_solver.solve(optical_workflow)
    for temperature in provider.electron_temperatures[previous_call_count:]:
        npt.assert_array_equal(temperature, [7200.0, 9900.0])
    npt.assert_array_equal(
        next_result.plasma_state.electron_temperature, [7200.0, 9900.0]
    )


@pytest.mark.parametrize(
    "include_heating, expected_change", [(False, 0.0), (True, 99.0)]
)
def test_nonthermal_heating_change_is_optional_for_convergence(
    include_heating: bool,
    expected_change: float,
) -> None:
    """A diagnostic heating change cannot block a fixed-temperature solve."""
    ionization = pd.DataFrame([[1.0e-5]], columns=[54])
    previous = FixedRateState(
        ionization_rates=ionization,
        heating_rate_density=pd.Series([1.0e-11], index=[54]),
    )
    candidate = FixedRateState(
        ionization_rates=ionization,
        heating_rate_density=pd.Series([1.0e-9], index=[54]),
    )
    change = NonthermalHeliumPlasmaSolver._rate_change(
        candidate, previous, include_heating=include_heating
    )
    npt.assert_allclose(change, expected_change, rtol=1e-14, atol=0.0)


@pytest.mark.parametrize(
    ("initial_population", "target_population", "accepted"),
    [
        (5.0e-21, 0.0, True),
        (2.0e-20, 0.0, False),
        (0.0, 5.0e-21, True),
        (0.0, 2.0e-20, False),
        (1.0e-20, 0.0, True),
        (1.0, 1.0005, True),
        (1.0005, 1.0, True),
        (1.0, 1.0010005, True),
        (1.0010005, 1.0, True),
        (1.0, 1.002, False),
        (1.0e-18, 1.0105e-18, True),
        (1.0e-18, 1.012e-18, False),
        (0.0, 0.0, True),
    ],
)
def test_nonthermal_level_acceptance_uses_absolute_and_relative_tolerances(
    initial_population: float,
    target_population: float,
    accepted: bool,
) -> None:
    """Verify mixed-tolerance acceptance independently of plasma physics.

    Regime: one trace level in cm^-3, otherwise fixed equilibrium and rates.
    The first solve preserves the seed; the second proposes the target.
    Verification: hand-selected differences bracket atol=1e-20 cm^-3 and
    rtol=1e-3. Small and populated limits are checked in both directions.
    This orchestration test detects treating atol as a denominator floor or
    accepting significant relative changes merely because a level is small.
    """

    class TraceLevelEquilibriumSolver(FixedEquilibriumSolver):
        def solve(
            self, *args: object, **kwargs: object
        ) -> PlasmaEquilibriumSolution:
            solution = super().solve(*args, **kwargs)
            populations = solution.evaluation.absolute_level_population
            populations.loc[(2, 0, 0)] = 6.0e8
            populations.loc[(2, 0, 1)] = (
                initial_population
                if len(self.level_seeds) == 1
                else target_population
            )
            return solution

    plasma = FakePlasma()
    plasma.electron_densities[:] = 4.0e8
    plasma.t_electrons[:] = 8000.0
    plasma.ion_number_density.loc[(2, 0)] = 6.0e8
    plasma.ion_number_density.loc[(2, 1)] = 4.0e8
    plasma.level_number_density.loc[(2, 0, 1)] = initial_population
    optical_workflow = SimpleNamespace(
        plasma_solver=plasma,
        simulation_state=SimpleNamespace(t_radiative=np.array([1.0e4])),
    )

    def evaluator_factory(*args: object) -> SimpleNamespace:
        del args
        return SimpleNamespace()

    coupling = NonthermalHeliumPlasmaSolver(
        FakeNonthermalSolver(),
        evaluator_factory,
        equilibrium_solver=TraceLevelEquilibriumSolver(),
        settings=NonthermalCouplingSettings(maximum_iterations=2),
    )
    coupling.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        pd.Index([0]),
        np.ones(1) * u.erg / u.s / u.cm**3,
        np.zeros(1) * u.erg / u.s / u.cm**3,
        np.ones(1) * u.erg / u.s / u.cm**3,
    )
    if accepted:
        coupling.solve(optical_workflow)
        npt.assert_allclose(
            plasma.level_number_density.loc[(2, 0, 1)],
            target_population,
            rtol=1e-14,
            atol=0.0,
        )
    else:
        with pytest.raises(RuntimeError, match="blocker=level_population"):
            coupling.solve(optical_workflow)
        assert coupling.last_result is None

    last = coupling.last_convergence_history.iloc[-1]
    assert last["failing_level_count"] == (0 if accepted else 1)
    npt.assert_allclose(
        coupling.last_iteration_diagnostics.loc["level_population", "change"],
        last["level_population"],
        rtol=1e-14,
        atol=0.0,
    )


@pytest.mark.parametrize("target_population", [0.0, 1.0e-30])
def test_nonthermal_selective_relaxation_removes_zero_target_tail(
    caplog: pytest.LogCaptureFixture,
    target_population: float,
) -> None:
    """Apply a vanishing level fraction without a geometric damping tail.

    Regime: 50%-relaxed charge/ion updates from a trace level at 1e-8 cm^-3
    to zero or 1e-30 cm^-3. Verification: only the first input has the large
    trace population; later inputs use the candidate level fraction scaled
    to the damped He I population. There is no twenty-step halving tail.
    This is an orchestration test; it does not validate helium atomic rates.
    """

    class ZeroTargetEquilibriumSolver(FixedEquilibriumSolver):
        def solve(
            self, *args: object, **kwargs: object
        ) -> PlasmaEquilibriumSolution:
            solution = super().solve(*args, **kwargs)
            populations = solution.evaluation.absolute_level_population
            populations.loc[(2, 0, 0)] = 6.0e8
            populations.loc[(2, 0, 1)] = target_population
            return solution

    plasma = FakePlasma()
    for populations in (
        plasma.number_density,
        plasma.ion_number_density,
        plasma.level_number_density,
    ):
        populations.columns = pd.Index([53])
    plasma.electron_densities.index = pd.Index([53])
    plasma.level_number_density.loc[(2, 0, 0)] = 8.0e8
    plasma.level_number_density.loc[(2, 0, 1)] = 1.0e-8
    optical_workflow = SimpleNamespace(
        plasma_solver=plasma,
        simulation_state=SimpleNamespace(t_radiative=np.array([1.0e4])),
    )

    def evaluator_factory(*args: object) -> SimpleNamespace:
        del args
        return SimpleNamespace()

    coupling_solver = NonthermalHeliumPlasmaSolver(
        FakeNonthermalSolver(),
        evaluator_factory,
        equilibrium_solver=ZeroTargetEquilibriumSolver(),
    )
    coupling_solver.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        pd.Index([53]),
        np.ones(1) * u.erg / u.s / u.cm**3,
        np.zeros(1) * u.erg / u.s / u.cm**3,
        np.ones(1) * u.erg / u.s / u.cm**3,
    )
    caplog.set_level(
        "INFO", logger="tardis.workflows.nonthermal_helium_workflow"
    )

    coupling_solver.solve(optical_workflow)

    history = coupling_solver.last_convergence_history
    assert len(history) < 20
    assert history.iloc[0]["zero_target_level_count"] == 1
    npt.assert_array_equal(history["zero_target_level_count"].iloc[1:], 0)
    assert history.iloc[-1]["electron_density"] < 1.0e-3
    assert history.iloc[-1]["level_population"] < 1.0e-3
    assert "level=(2, 0, 1)" in caplog.text
    assert "shell=53" in caplog.text
    assert "zero/below-atol targets=1" in caplog.text
    assert coupling_solver.last_zero_target_levels.empty
    assert coupling_solver.last_result is not None
    npt.assert_allclose(
        plasma.level_number_density.loc[(2, 0, 1)],
        target_population,
        rtol=1e-14,
        atol=0.0,
    )


@pytest.mark.parametrize(
    ("initial_he_ii_density", "target_he_ii_density"),
    [(2.0e8, 4.0e8), (2.0e8, 0.0), (0.0, 0.0)],
)
def test_nonthermal_selective_relaxation_conserves_ion_populations(
    initial_he_ii_density: float,
    target_he_ii_density: float,
) -> None:
    """Preserve ion totals while replacing level fractions without damping.

    Regime: two helium ions, two levels per ion, two labelled shells, alpha=0.5.
    Verification: the next input has midpoint ion/electron/temperature values,
    candidate level fractions, and level sums equal to each damped ion total.
    For a vanished candidate ion, its surviving damped population retains the
    previous shell-specific fractions. Exact analytic fractions are compared
    at rtol=1e-14, allowing only floating-point roundoff. This verifies the
    update/conservation law, not the physical rate solution.
    """
    plasma = FakePlasma(shell_count=2)
    shell_columns = pd.Index([25, 53])
    plasma.number_density.columns = shell_columns
    plasma.ion_number_density.columns = shell_columns
    plasma.electron_densities = pd.Series(
        [initial_he_ii_density] * 2, index=shell_columns
    )
    plasma.ion_number_density.loc[(2, 0)] = 1.0e9 - initial_he_ii_density
    plasma.ion_number_density.loc[(2, 1)] = initial_he_ii_density
    plasma.t_electrons[:] = 6000.0
    level_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0), (2, 0, 1), (2, 1, 0), (2, 1, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    plasma.level_number_density = pd.DataFrame(
        [
            [0.75 * (1.0e9 - initial_he_ii_density)] * 2,
            [0.25 * (1.0e9 - initial_he_ii_density)] * 2,
            [0.75 * initial_he_ii_density, 0.25 * initial_he_ii_density],
            [0.25 * initial_he_ii_density, 0.75 * initial_he_ii_density],
        ],
        index=level_index,
        columns=shell_columns,
    )
    original_levels = plasma.level_number_density.copy(deep=True)
    captured_states = []

    class TwoIonEquilibriumSolver(FixedEquilibriumSolver):
        def solve(
            self, *args: object, **kwargs: object
        ) -> PlasmaEquilibriumSolution:
            solution = super().solve(*args, **kwargs)
            solution.evaluation.ion_population.loc[(2, 0)] = (
                1.0e9 - target_he_ii_density
            )
            solution.evaluation.ion_population.loc[(2, 1)] = (
                target_he_ii_density
            )
            solution.evaluation.absolute_level_population = pd.DataFrame(
                np.tile(
                    [
                        [0.25 * (1.0e9 - target_he_ii_density)],
                        [0.75 * (1.0e9 - target_he_ii_density)],
                        [target_he_ii_density],
                        [0.0],
                    ],
                    (1, 2),
                ),
                index=level_index,
                columns=shell_columns,
            )
            return PlasmaEquilibriumSolution(
                np.full(2, target_he_ii_density),
                solution.electron_temperature,
                solution.evaluation,
            )

    def evaluator_factory(
        state: CoupledPlasmaState, *args: object
    ) -> SimpleNamespace:
        del args
        captured_states.append(state)
        return SimpleNamespace()

    coupling = NonthermalHeliumPlasmaSolver(
        FakeNonthermalSolver(),
        evaluator_factory,
        equilibrium_solver=TwoIonEquilibriumSolver(),
        settings=NonthermalCouplingSettings(maximum_iterations=2),
    )
    coupling.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        shell_columns,
        np.ones(2) * u.erg / u.s / u.cm**3,
        np.zeros(2) * u.erg / u.s / u.cm**3,
        np.ones(2) * u.erg / u.s / u.cm**3,
    )
    with pytest.raises(RuntimeError, match="did not converge"):
        coupling.solve(
            SimpleNamespace(
                plasma_solver=plasma,
                simulation_state=SimpleNamespace(t_radiative=np.full(2, 1.0e4)),
            )
        )
    next_state = captured_states[1]
    expected_he_i = 1.0e9 - 0.5 * (initial_he_ii_density + target_he_ii_density)
    expected_he_ii = 0.5 * (initial_he_ii_density + target_he_ii_density)
    npt.assert_allclose(next_state.electron_density, expected_he_ii, rtol=1e-14)
    npt.assert_allclose(next_state.electron_temperature, 7000.0, rtol=1e-14)
    npt.assert_allclose(
        next_state.ion_number_density,
        [[expected_he_i] * 2, [expected_he_ii] * 2],
        rtol=1e-14,
    )
    pdt.assert_frame_equal(
        next_state.level_number_density.groupby(
            level=["atomic_number", "ion_number"]
        ).sum(),
        next_state.ion_number_density,
        rtol=1e-14,
    )
    npt.assert_allclose(
        next_state.level_number_density.loc[(2, 0)] / expected_he_i,
        [[0.25, 0.25], [0.75, 0.75]],
        rtol=1e-14,
    )
    if expected_he_ii == 0.0:
        npt.assert_array_equal(next_state.level_number_density.loc[(2, 1)], 0.0)
    else:
        npt.assert_allclose(
            next_state.level_number_density.loc[(2, 1)] / expected_he_ii,
            (
                [[1.0, 1.0], [0.0, 0.0]]
                if target_he_ii_density > 0.0
                else [[0.75, 0.25], [0.25, 0.75]]
            ),
            rtol=1e-14,
            atol=0.0,
        )
    assert coupling.last_result is None
    pdt.assert_frame_equal(plasma.level_number_density, original_levels)


def test_nonthermal_solver_uses_workflow_evaluator_factory_by_default() -> None:
    """Keep continuum-estimator ownership inside the optical workflow.

    Claim
    -----
    A caller does not need to reconstruct continuum mappings in a notebook.

    Regime
    ------
    One-shell deterministic coupling with no explicit evaluator factory.

    Verification
    ------------
    The coupling solver requests every evaluator from the optical workflow.
    """
    provider = FakeNonthermalSolver()
    equilibrium_solver = FixedEquilibriumSolver()
    factory_calls = []
    optical_workflow = SimpleNamespace(
        plasma_solver=FakePlasma(),
        simulation_state=SimpleNamespace(t_radiative=np.array([1.0e4])),
    )

    def build_helium_evaluator(
        state: object,
        fixed_rates: object,
        bound_bound_species_filter: tuple[int, int] | None,
    ) -> SimpleNamespace:
        factory_calls.append((state, fixed_rates, bound_bound_species_filter))
        return SimpleNamespace()

    optical_workflow.build_helium_evaluator = build_helium_evaluator
    coupling_solver = NonthermalHeliumPlasmaSolver(
        provider,
        evaluator_factory=None,
        equilibrium_solver=equilibrium_solver,
        settings=NonthermalCouplingSettings(
            maximum_iterations=4,
            relative_tolerance=1.0e-12,
            relaxation=1.0,
        ),
    )
    coupling_solver.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        pd.Index([0]),
        np.array([1.0]) * u.erg / u.s / u.cm**3,
        np.array([0.0]) * u.erg / u.s / u.cm**3,
        np.array([1.0]) * u.erg / u.s / u.cm**3,
    )

    coupling_solver.solve(optical_workflow)

    assert len(factory_calls) == 3
    assert all(call[2] == (2, 0) for call in factory_calls)


def test_nonthermal_solver_updates_all_shells() -> None:
    """The coupled solve evaluates and commits every simulation shell.

    Claim
    -----
    Every shell participates in the coupled nonthermal equilibrium solve.

    Regime
    ------
    A deterministic 55-shell plasma spans the former development cutoff.

    Verification
    ------------
    The evaluator receives all physical shell labels and every shell equals
    the exact equilibrium state returned by the test solver.
    """
    provider = FakeNonthermalSolver()
    equilibrium_solver = FixedEquilibriumSolver()
    plasma = FakePlasma(shell_count=55)
    factory_shell_columns = []
    optical_workflow = SimpleNamespace(
        plasma_solver=plasma,
        simulation_state=SimpleNamespace(
            t_radiative=np.linspace(1.0e4, 1.1e4, 55)
        ),
    )

    def evaluator_factory(
        state: CoupledPlasmaState,
        fixed_rates: object,
        workflow: object,
        bound_bound_species_filter: tuple[int, int] | None,
    ) -> SimpleNamespace:
        del fixed_rates, workflow, bound_bound_species_filter
        factory_shell_columns.append(state.ion_number_density.columns.copy())
        return SimpleNamespace()

    coupling_solver = NonthermalHeliumPlasmaSolver(
        provider,
        evaluator_factory,
        equilibrium_solver=equilibrium_solver,
        settings=NonthermalCouplingSettings(
            maximum_iterations=4,
            relative_tolerance=1.0e-12,
            relaxation=1.0,
        ),
    )
    coupling_solver.deposition_state = DepositionState(
        1.0 * u.s,
        0,
        pd.RangeIndex(55),
        np.ones(55) * u.erg / u.s / u.cm**3,
        np.zeros(55) * u.erg / u.s / u.cm**3,
        np.ones(55) * u.erg / u.s / u.cm**3,
    )

    result = coupling_solver.solve(optical_workflow)

    assert all(
        columns.equals(pd.RangeIndex(55)) for columns in factory_shell_columns
    )
    npt.assert_array_equal(
        result.plasma_state.electron_density, np.full(55, 4.0e8)
    )
    npt.assert_array_equal(
        result.plasma_state.electron_temperature,
        np.full(55, 8000.0),
    )
    pdt.assert_frame_equal(
        result.plasma_state.ion_number_density,
        pd.DataFrame(
            np.tile([[6.0e8], [4.0e8]], 55),
            index=result.plasma_state.ion_number_density.index,
            columns=pd.RangeIndex(55),
        ),
    )
    pdt.assert_frame_equal(
        result.plasma_state.level_number_density,
        pd.DataFrame(
            np.tile([[4.0e8], [2.0e8]], 55),
            index=result.plasma_state.level_number_density.index,
            columns=pd.RangeIndex(55),
        ),
    )


class FakeHighEnergyWorkflow:
    """Return a deterministic flight-path deposition result."""

    def __init__(self, outer_velocity: float = 2.0e8) -> None:
        self.simulation_state = SimpleNamespace(
            v_inner=np.array([1.0e8]) * u.cm / u.s,
            v_outer=np.array([outer_velocity]) * u.cm / u.s,
        )
        self.run_arguments: dict[str, object] | None = None

    def run(self, **kwargs: object) -> SimpleNamespace:
        self.run_arguments = kwargs
        return SimpleNamespace(
            gamma_ray_deposition_estimator=pd.DataFrame(
                [[2.0]], index=[0], columns=[86400.0]
            ),
            positron_energy=pd.DataFrame([[0.0]], index=[0], columns=[86400.0]),
        )

    def get_times(
        self,
        time_start: float,
        time_end: float,
        time_space: str,
        time_steps: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        del time_start, time_end, time_space, time_steps
        return np.array([1.0, 2.0]), np.array([1.5])


class FakeOpticalWorkflow:
    """Record execution after a coupling solver has been attached."""

    def __init__(self) -> None:
        self.plasma_solver = FakePlasma()
        self.simulation_state = SimpleNamespace(
            t_radiative=np.array([1.0e4]),
            time_explosion=1.5 * u.d,
            v_inner=np.array([1.0e8]) * u.cm / u.s,
            v_outer=np.array([2.0e8]) * u.cm / u.s,
        )
        self.ran = False

    def run(self) -> None:
        assert self.plasma_coupling_solver.deposition_state is not None
        self.ran = True


def make_high_energy_settings() -> HighEnergySettings:
    """Return the smallest complete high-energy argument set."""
    return HighEnergySettings(
        time_start=1.0,
        time_end=2.0,
        number_of_packets=1,
        time_steps=1,
        time_space="linear",
        seed=1963,
        positronium_fraction=0.0,
        spectrum_bins=1,
    )


def test_composed_workflow_runs_high_energy_once_and_attaches_coupling() -> (
    None
):
    """Prepare aligned deposition before entering optical transport."""
    optical_workflow = FakeOpticalWorkflow()
    high_energy_workflow = FakeHighEnergyWorkflow()
    coupling_solver = SimpleNamespace(deposition_state=None)
    workflow = NonthermalHeliumWorkflow(
        optical_workflow,
        high_energy_workflow,
        coupling_solver,
        make_high_energy_settings(),
    )

    workflow.run()

    assert optical_workflow.ran
    assert high_energy_workflow.run_arguments["seed"] == 1963
    npt.assert_allclose(
        coupling_solver.deposition_state.rate_density.to_value(
            u.erg / u.s / u.cm**3
        ),
        [2.0],
    )


def test_composed_workflow_rejects_mismatched_velocity_grid() -> None:
    """Do not transfer deposition between physically different shells."""
    workflow = NonthermalHeliumWorkflow(
        FakeOpticalWorkflow(),
        FakeHighEnergyWorkflow(outer_velocity=3.0e8),
        SimpleNamespace(deposition_state=None),
        make_high_energy_settings(),
    )

    with np.testing.assert_raises_regex(ValueError, "v_outer"):
        workflow.run()


def test_specialized_optical_workflow_uses_continuum_macro_atom() -> None:
    """Pass accepted continuum rates into macro-atom opacity construction.

    Claim
    -----
    The next packet iteration uses the accepted continuum state.

    Regime
    ------
    One deterministic opacity update with a continuum macro-atom solver.

    Verification
    ------------
    The macro-atom solver receives continuum rates and electron density.
    """
    macro_atom_arguments = []
    workflow = object.__new__(NonthermalHeliumOpticalWorkflow)
    workflow.plasma_solver = SimpleNamespace(
        j_blues=pd.DataFrame([[1.0]]),
        stimulated_emission_factor=np.ones((1, 1)),
        electron_densities=pd.Series([1.0]),
        t_electrons=np.array([1.0]),
    )
    workflow.atom_data = SimpleNamespace(
        lines=pd.DataFrame(
            {"nu": [1.0]},
            index=pd.Index([0]),
        )
    )
    continuum_opacity = object()
    continuum_rates = object()
    opacity_arguments = []
    workflow.continuum_opacity_state = continuum_opacity
    workflow.continuum_macro_atom_state = continuum_rates
    workflow.opacity_solver = SimpleNamespace(
        disable_line_scattering=False,
        solve=lambda *args, **kwargs: (
            opacity_arguments.append((args, kwargs)),
            SimpleNamespace(beta_sobolev=pd.DataFrame([[1.0]])),
        )[1],
    )
    workflow._update_sobolev_state = lambda: None
    workflow._tau_sobolev = pd.DataFrame([[1.0]])
    workflow._beta_sobolev = pd.DataFrame([[1.0]])
    workflow.macro_atom_solver = SimpleNamespace(
        solve=lambda *args: (macro_atom_arguments.append(args), object())[1]
    )

    workflow.solve_opacity()

    # Continuum opacity is now a separate transport state, not a line input.
    assert opacity_arguments[0][0] == (workflow.plasma_solver,)
    assert macro_atom_arguments[0][3] is continuum_rates
    assert (
        macro_atom_arguments[0][4] is workflow.plasma_solver.electron_densities
    )


@pytest.fixture
def helium_optical_workflow(
    tardis_regression_path: Path,
) -> NonthermalHeliumOpticalWorkflow:
    """Construct the three-shell H/He optical owner without Monte Carlo."""
    config = Configuration.from_yaml(
        "tardis/workflows/tests/data/iip_hydrogen_helium.yml"
    )
    config.atom_data = (
        tardis_regression_path
        / "atom_data"
        / "christians_atomdata_converted_04Dec25.h5"
    )
    # The classic bootstrap has no continuum; the specialized owner adds He I.
    config.plasma.continuum_interaction.species = []
    config.montecarlo.nthreads = 1
    config.plasma.continuum_interaction.nthreads = 2
    return NonthermalHeliumOpticalWorkflow(
        config,
        continuum_interaction_species=("He I",),
        show_progress_bars=False,
        log_level="ERROR",
    )


def test_helium_workflow_uses_explicit_continuum_mode(
    helium_optical_workflow: NonthermalHeliumOpticalWorkflow,
) -> None:
    """Keep continuum traversal explicit despite the global full-mode default.

    Orchestration claim: the actual workflow factory chooses the upstream
    explicit continuum path and independent opacity-worker limit, not optical
    Monte Carlo threads. The fixture specifies two opacity workers and one
    transport worker; exact comparisons detect default-full regressions.
    """
    workflow = helium_optical_workflow
    solver = workflow.macro_atom_solver
    assert (
        workflow.configuration.plasma.continuum_interaction.markov_chain.mode
        == "full"
    )
    assert solver.markov_chain_mode == "explicit"
    assert solver.markov_chain_elements == ()
    assert not solver.markov_chain_include_thermal_pool
    assert not solver.markov_chain_include_ionization_pool
    assert solver.nthreads == 2
    assert workflow.configuration.montecarlo.nthreads == 1


def test_helium_explicit_transport_refreshes_evaluator_estimators(
    helium_optical_workflow: NonthermalHeliumOpticalWorkflow,
) -> None:
    """Carry two actual transport steps into fixed-estimator helium rates.

    Orchestration verification on the three-shell real-data H/He fixture,
    without the high-energy workflow or a nonlinear plasma convergence run.
    Check the Lucy path-length photoionization normalization independently:
    raw estimator / (h t V), multiplied by the target-intensity damping.
    rtol=1e-12 permits only arithmetic roundoff in that normalization. This
    catches stale transport estimators and regression to reduced traversal;
    it is not an end-to-end nonthermal-population equivalence test.
    """
    workflow = helium_optical_workflow
    previous_estimators = None
    for iteration in range(2):
        workflow.completed_iterations = iteration
        opacity_states = workflow.solve_opacity()
        macro_state = opacity_states["macro_atom_state"]
        assert macro_state.absorbing_probability_matrix.shape == (3, 0, 0)
        npt.assert_array_equal(macro_state.macro_level_to_reduced_idx, -1)
        workflow.solve_montecarlo(opacity_states, no_of_real_packets=20000)
        workflow._capture_continuum_estimators()
        assert workflow._continuum_estimators is not previous_estimators
        previous_estimators = workflow._continuum_estimators

        state = workflow._current_plasma_state()
        evaluator = workflow.build_helium_evaluator(
            state, FixedRateState(), (2, 0)
        )
        rate_solver = evaluator.ion_population_solver.rate_matrix_solver.radiative_ionization_rate_solver
        actual, _, _ = rate_solver.solve_coefficients(
            state.electron_temperature * u.K
        )
        transport = workflow.transport_state
        volume = transport.geometry_state_numba.volume
        time = transport.time_of_simulation.value
        target_intensity = (
            workflow.simulation_state.dilution_factor
            * workflow.simulation_state.t_radiative.value**4
            * const.sigma_sb.cgs.value
            / np.pi
        )
        measured_intensity = transport.estimators_bulk.mean_intensity_total / (
            4.0 * np.pi * time * volume
        )
        damping = target_intensity / measured_intensity
        expected = (
            transport.estimators_continuum.photo_ion_estimator
            / (const.h.cgs.value * time * volume)
            * damping
        )
        # Continuum arrays are reordered by their declared edge identities.
        edge_positions = workflow.plasma_solver.level2continuum_idx.to_numpy(
            dtype=int
        )
        assert np.any(expected > 0.0)
        npt.assert_allclose(
            np.asarray(actual), expected[edge_positions], rtol=1e-12, atol=0.0
        )


def test_helium_evaluator_constructs_nebular_background(
    helium_optical_workflow: NonthermalHeliumOpticalWorkflow,
) -> None:
    """Build a real helium evaluator without running high-energy or optical MC.

    This is orchestration verification on the three-shell H/He fixture. The
    hydrogen background must receive the radiation-field nebular factor,
    independently evaluated through the initial plasma graph, rather than
    thermal LTE. Dilute excitation tables and physical shell labels must also
    reach the evaluator. rtol=1e-12 allows roundoff in identical physical inputs.
    """
    workflow = helium_optical_workflow
    state = workflow._current_plasma_state()
    evaluator = workflow.build_helium_evaluator(state, FixedRateState(), (2, 0))
    background = evaluator.dilute_radiation
    assert background is not None
    npt.assert_array_equal(
        background.radiation_temperature,
        np.asarray(workflow.plasma_solver.t_rad),
    )
    npt.assert_array_equal(
        background.dilution_factor, np.asarray(workflow.plasma_solver.w)
    )
    pdt.assert_frame_equal(
        background.level_boltzmann_factor,
        workflow.plasma_solver.level_boltzmann_factor,
    )
    nebular = evaluator.calculate_nebular_ionization_factor(
        state.electron_temperature
    )
    npt.assert_allclose(
        nebular.loc[(1, 1)],
        workflow.plasma_solver.phi.loc[(1, 1)],
        rtol=1e-12,
        atol=0.0,
    )


def test_helium_evaluator_initializes_current_adjacent_ion_budget(
    helium_optical_workflow: NonthermalHeliumOpticalWorkflow,
) -> None:
    """Exclude prescribed He III from the reduced adjacent-ion density budget.

    Claim: current He I + He II, not stale elemental or ion totals, sets the
    reduced geometry. Regime: conserved three-stage helium in three shells,
    including a He III-free limit and reversed shell order. Verification uses
    N_He - prescribed He III independently of the pair sum; rtol=1e-14 allows
    addition/subtraction roundoff only.
    """
    workflow = helium_optical_workflow
    original = workflow._current_plasma_state()
    total_helium = workflow.plasma_solver.number_density.loc[2]
    for shell_order in (
        original.ion_number_density.columns,
        total_helium.index[::-1],
    ):
        for he_iii_fraction in (
            np.array([0.1, 0.25, 0.4]),
            np.array([0.3, 0.45, 0.6]),
            np.zeros(3),
        ):
            current_ions = original.ion_number_density.copy(deep=True)
            prescribed_he_iii = total_helium * he_iii_fraction
            current_ions.loc[(2, 0)] = total_helium * 0.2
            current_ions.loc[(2, 1)] = total_helium * (0.8 - he_iii_fraction)
            current_ions.loc[(2, 2)] = prescribed_he_iii
            state = CoupledPlasmaState(
                original.electron_density,
                original.electron_temperature,
                current_ions.loc[:, shell_order],
                original.level_number_density.loc[:, shell_order],
            )
            evaluator = workflow.build_helium_evaluator(
                state, FixedRateState(), (2, 0)
            )
            npt.assert_allclose(
                [
                    geometry.species_number_density
                    for geometry in evaluator.shell_number_densities
                ],
                (total_helium - prescribed_he_iii).loc[shell_order],
                rtol=1e-14,
                atol=0.0,
                err_msg="Reduced He I/He II geometry includes an unrelated ion stage.",
            )


def test_helium_evaluator_preserves_current_population_geometry(
    helium_optical_workflow: NonthermalHeliumOpticalWorkflow,
) -> None:
    """Preserve caller-owned populations and level ordering during construction."""
    workflow = helium_optical_workflow
    state = workflow._current_plasma_state()
    original_ions = state.ion_number_density.copy(deep=True)
    original_levels = state.level_number_density.copy(deep=True)
    evaluator = workflow.build_helium_evaluator(state, FixedRateState(), (2, 0))

    pdt.assert_frame_equal(state.ion_number_density, original_ions)
    pdt.assert_frame_equal(state.level_number_density, original_levels)
    pdt.assert_frame_equal(
        workflow.plasma_solver.ion_number_density, original_ions
    )
    pdt.assert_frame_equal(
        workflow.plasma_solver.level_number_density, original_levels
    )
    expected_positions = np.flatnonzero(
        (original_levels.index.get_level_values("atomic_number") == 2)
        & (original_levels.index.get_level_values("ion_number") == 0)
    )
    for shell, geometry in zip(
        original_levels.columns, evaluator.shell_number_densities, strict=True
    ):
        npt.assert_array_equal(
            geometry.level_number_density, original_levels[shell]
        )
        npt.assert_array_equal(
            geometry.species_level_positions, expected_positions
        )


def test_specialized_optical_workflow_enables_continuum_transport() -> None:
    """Provide the IIP transport with continuum opacity and macro-atom data."""
    initialized_with = []
    transport_state = SimpleNamespace(
        packet_collection=SimpleNamespace(output_energies=np.array([1.0]))
    )
    workflow = object.__new__(NonthermalHeliumOpticalWorkflow)
    workflow.simulation_state = object()
    workflow.continuum_opacity_state = object()
    workflow.completed_iterations = 0
    workflow.show_progress_bars = False
    workflow.transport_solver = SimpleNamespace(
        initialize_transport_state=lambda *args, **kwargs: (
            initialized_with.append(args),
            transport_state,
        )[1],
        run=lambda *args, **kwargs: np.array([2.0]),
    )

    opacity_state = object()
    macro_atom_state = object()
    virtual_packet_energies = workflow.solve_montecarlo(
        {"opacity_state": opacity_state, "macro_atom_state": macro_atom_state},
        no_of_real_packets=1,
    )

    assert initialized_with[0][1] is opacity_state
    assert initialized_with[0][2] is macro_atom_state
    assert initialized_with[0][3] is workflow.continuum_opacity_state
    assert initialized_with[0][4] == 1
    assert virtual_packet_energies is None


@pytest.mark.parametrize("wavelength_angstrom", [584.0, 5000.0, 20000.0])
@pytest.mark.parametrize("line_estimator", [0.0, 2.0])
@pytest.mark.parametrize("intensity_damping", [0.5, 2.0])
def test_line_intensity_damps_packets_but_not_target_fallback(
    wavelength_angstrom: float,
    line_estimator: float,
    intensity_damping: float,
) -> None:
    """Verify measured J and target Planck/epsilon branches independently.

    The optical window uses c t_exp/(4 pi t_sim V) times the packet sum
    and its bolometric damping. UV/IR substitutions and zero-packet epsilon
    floors already use the final target's w B_nu and must not be damped.
    Analytic normalization and Planck's law provide independent oracles;
    rtol=1e-13 permits floating-point arithmetic, with no absolute floor.
    """
    temperature = 8000.0
    dilution_factor = 0.3
    frequency = (wavelength_angstrom * u.AA).to_value(u.Hz, u.spectral())
    target_intensity = (
        dilution_factor * const.sigma_sb.cgs.value * temperature**4 / np.pi
    )
    time_simulation = 2.0 * u.s
    time_explosion = 3.0 * u.s
    volume = np.array([5.0])
    workflow = object.__new__(NonthermalHeliumOpticalWorkflow)
    workflow.transport_state = SimpleNamespace(
        estimators_continuum=init_estimators_continuum((1, 1), 1),
        geometry_state_numba=SimpleNamespace(volume=volume),
        time_of_simulation=time_simulation,
        time_explosion=time_explosion,
        estimators_bulk=SimpleNamespace(
            mean_intensity_total=(
                target_intensity
                * 4.0
                * np.pi
                * time_simulation.value
                * volume
                / intensity_damping
            )
        ),
        estimators_line=SimpleNamespace(
            mean_intensity_blueward=np.array([[line_estimator]])
        ),
        opacity_state_numba=SimpleNamespace(line_list_nu=np.array([frequency])),
    )
    radiation_field = DilutePlanckianRadiationField(
        np.array([temperature]) * u.K, np.array([dilution_factor])
    )
    workflow.simulation_state = SimpleNamespace(
        dilution_factor=np.array([dilution_factor]),
        t_radiative=np.array([temperature]) * u.K,
        radiation_field_state=radiation_field,
    )
    updated = []
    workflow.plasma_solver = SimpleNamespace(
        number_density=pd.DataFrame([[1.0]]),
        ion_number_density=pd.DataFrame(
            [[1.0]],
            index=pd.MultiIndex.from_tuples(
                [(2, 1)], names=["atomic_number", "ion_number"]
            ),
        ),
        electron_densities=pd.Series([1.0]),
        t_electrons=np.array([temperature]),
        lines=pd.DataFrame(index=pd.Index([0])),
        update=lambda **kwargs: updated.append(kwargs),
    )
    radiation_solver = MCRadiationFieldPropertiesSolver(w_epsilon=1.0e-10)
    workflow.transport_solver = SimpleNamespace(
        radfield_prop_solver=radiation_solver
    )

    workflow._capture_continuum_estimators()

    planck_intensity = (
        dilution_factor
        * 2.0
        * const.h.cgs.value
        * frequency**3
        / const.c.cgs.value**2
        / np.expm1(
            const.h.cgs.value * frequency / (const.k_B.cgs.value * temperature)
        )
    )
    if line_estimator == 0.0:
        expected = radiation_solver.w_epsilon * planck_intensity
    elif 2500.0 < wavelength_angstrom < 10000.0:
        expected = (
            line_estimator
            * const.c.cgs.value
            * time_explosion.value
            / (4.0 * np.pi * time_simulation.value * volume[0])
            * intensity_damping
        )
    else:
        expected = planck_intensity
    npt.assert_allclose(
        updated[0]["j_blues"].iloc[0, 0], expected, rtol=1.0e-13, atol=0.0
    )


def test_continuum_estimators_remain_finite_in_unvisited_shells() -> None:
    """An unvisited shell keeps zero packet rates without NaN damping.

    Claim
    -----
    Sparse Monte Carlo estimators do not contaminate the plasma rates.

    Regime
    ------
    Two shells with radiation-field paths only in the first shell.

    Verification
    ------------
    Every normalized continuum estimator and line intensity is finite, and
    the unvisited shell's continuum rate remains exactly zero.
    """
    workflow = object.__new__(NonthermalHeliumOpticalWorkflow)
    shell_columns = pd.RangeIndex(2)
    raw = init_estimators_continuum((1, 2), 2)
    raw.photo_ion_estimator[0, 0] = 2.0
    workflow.transport_state = SimpleNamespace(
        estimators_continuum=raw,
        geometry_state_numba=SimpleNamespace(volume=np.ones(2)),
        time_of_simulation=1.0 * u.s,
        time_explosion=1.0 * u.s,
        estimators_bulk=SimpleNamespace(
            mean_intensity_total=np.array([1.0, 0.0])
        ),
        estimators_line=SimpleNamespace(
            mean_intensity_blueward=np.zeros((1, 2))
        ),
        opacity_state_numba=SimpleNamespace(line_list_nu=np.array([1.0])),
    )
    workflow.simulation_state = SimpleNamespace(
        dilution_factor=np.ones(2),
        t_radiative=np.array([6000.0, 6000.0]) * u.K,
        radiation_field_state=object(),
    )
    updated = []
    workflow.plasma_solver = SimpleNamespace(
        number_density=pd.DataFrame({0: [1.0], 1: [1.0]}),
        ion_number_density=pd.DataFrame(
            [[1.0, 1.0]],
            index=pd.MultiIndex.from_tuples(
                [(2, 1)], names=["atomic_number", "ion_number"]
            ),
            columns=shell_columns,
        ),
        electron_densities=pd.Series([1.0, 1.0]),
        t_electrons=np.array([1.0, 1.0]),
        lines=pd.DataFrame(index=pd.Index([0])),
        update=lambda **kwargs: updated.append(kwargs),
    )
    workflow.transport_solver = SimpleNamespace(
        radfield_prop_solver=SimpleNamespace(
            estimate_jblues=lambda *args, **kwargs: np.array([[1.0, 0.0]])
        )
    )

    workflow._capture_continuum_estimators()

    normalized = workflow._continuum_estimators
    assert np.isfinite(normalized.photo_ion_estimator).all()
    assert np.isfinite(normalized.ff_heating_estimator).all()
    assert np.isfinite(updated[0]["j_blues"].to_numpy()).all()
    assert normalized.photo_ion_estimator[0, 1] == 0.0


def test_specialized_optical_workflow_couples_after_estimator_capture(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Build next continuum state after accepting coupled helium plasma."""
    call_order = []
    workflow = object.__new__(NonthermalHeliumOpticalWorkflow)
    estimated_radfield_properties = object()
    monkeypatch.setattr(
        StandardTARDISWorkflow,
        "solve_plasma",
        lambda owner, estimated: call_order.append(
            ("classic", owner, estimated)
        ),
    )
    workflow._update_sobolev_state = lambda: call_order.append(("sobolev",))
    workflow._capture_continuum_estimators = lambda: call_order.append(
        ("estimators",)
    )
    workflow._build_continuum_states = lambda *args: call_order.append(
        ("continuum", *args)
    )
    coefficients = object()
    saha_factor = object()
    workflow.plasma_coupling_solver = SimpleNamespace(
        solve=lambda owner: (
            call_order.append(("nonthermal", owner)),
            SimpleNamespace(
                equilibrium_solution=SimpleNamespace(
                    evaluation=SimpleNamespace(
                        continuum_coefficients=coefficients,
                        level_to_continuum_saha_factor=saha_factor,
                    )
                )
            ),
        )[1]
    )

    workflow.solve_plasma(estimated_radfield_properties)

    assert [call[0] for call in call_order] == [
        "classic",
        "estimators",
        "nonthermal",
        "continuum",
        "sobolev",
    ]
    assert call_order[0][2] is estimated_radfield_properties
    assert call_order[3][1:] == (coefficients, saha_factor)
