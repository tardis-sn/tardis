"""Distinguish open inner kinetics from changes in the terminal state."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.equilibrium.diagnostics import diagnose_level_closure
from tardis.plasma.equilibrium.evaluator import PlasmaEquilibriumEvaluator
from tardis.plasma.equilibrium.inputs import (
    BoundBoundMatrixRates,
    ContinuumRateCoefficients,
    ShellNumberDensity,
)
from tardis.plasma.equilibrium.tests.test_evaluator import (
    toy_evaluator as toy_evaluator,  # noqa: PLC0414
)


@pytest.fixture
def closure_context(
    toy_evaluator: PlasmaEquilibriumEvaluator, monkeypatch: pytest.MonkeyPatch
) -> SimpleNamespace:
    """Supply analytic rates and a deliberately prescribed ionic boundary.

    At ne=1, gamma=[1,1], alpha=[1,1] and downward rate2 give
    f=[5/6,1/6], eliminated upper/lower ratio1/2. Prescribing ions2/2
    instead isolates an inconsistent final ratio without a global ion solve.
    This is not a claimed self-consistent charge solution.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()
    evaluator.elemental_number_density.columns = pd.Index([17])
    evaluator.elemental_number_density.iloc[0, 0] = 4.0
    evaluator.j_blues.columns = pd.Index([17])
    evaluator.shell_number_densities = (
        ShellNumberDensity(4.0, np.array([4.0, 0.0]), np.array([0, 1])),
    )
    evaluator.sobolev_inputs = (
        replace(evaluator.sobolev_inputs[0], tau_coefficient=np.zeros(1)),
    )
    context = SimpleNamespace(
        evaluator=evaluator,
        fractions=np.array([5 / 6, 1 / 6]),
        final_density=1.0,
        ions=pd.DataFrame(
            [[2.0], [2.0]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
            ),
            columns=[17],
        ),
        bound=BoundBoundMatrixRates(
            2,
            np.array([1]),
            np.array([0]),
            np.zeros(1),
            np.zeros(1),
            np.array([2.0]),
            np.array([-1]),
        ),
        continuum=ContinuumRateCoefficients(
            np.ones(2), np.zeros(2), np.ones(2), np.zeros(2), np.zeros(2)
        ),
    )
    original_coefficients = evaluator.calculate_continuum_coefficients

    def continuum_coefficients(temperature: npt.NDArray[np.float64]) -> tuple:
        coefficients = original_coefficients(temperature)
        return ((context.continuum,), *coefficients[1:])

    def level_solution(
        shell_idx: int,
        density: float,
        temperature: float,
        seed: npt.NDArray[np.float64],
        continuum: ContinuumRateCoefficients,
        geometry: ShellNumberDensity,
        sobolev: object,
        bound: BoundBoundMatrixRates,
    ) -> tuple:
        residual, beta, ratio = evaluator._calculate_level_state(
            shell_idx,
            density,
            temperature,
            context.fractions,
            continuum,
            geometry,
            sobolev,
            bound_bound_rates=bound,
        )
        return context.fractions.copy(), ratio, beta, residual

    monkeypatch.setattr(
        evaluator, "calculate_continuum_coefficients", continuum_coefficients
    )
    monkeypatch.setattr(
        evaluator.rate_matrix_solver,
        "prepare_species_rates",
        lambda *args: (context.bound,),
    )
    monkeypatch.setattr(evaluator, "_calculate_level_solution", level_solution)
    monkeypatch.setattr(
        evaluator,
        "_solve_charge_population",
        lambda *args: (
            context.ions.copy(deep=True),
            pd.Series([context.final_density], index=[17]),
        ),
    )
    return context


def test_baseline_closed_inner_can_have_open_authoritative_ratio(
    closure_context: SimpleNamespace,
) -> None:
    """Existing equations close with r=.5 but gain source .5/level at r=1.

    This pre-diagnostic verification establishes the physical distinction;
    a missing new API is not being counted as scientific red-phase evidence.
    """
    context = closure_context
    result = context.evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    np.testing.assert_allclose(
        result.trial_level_residual, 0.0, atol=2e-16, rtol=0.0
    )
    np.testing.assert_allclose(
        result.diagnostic_ion_ratio, 0.5, atol=0.0, rtol=1e-15
    )
    allowance = 64 * np.finfo(float).eps / 1e-8 * (11 / 3)
    np.testing.assert_allclose(
        result.relative_level_residual.to_numpy().ravel(),
        [0.5 / (13 / 6 + allowance), 0.5 / (3 / 2 + allowance)],
        atol=0.0,
        rtol=1e-13,
    )


@pytest.mark.parametrize("inner_open", [False, True])
def test_diagnostic_distinguishes_inner_and_ion_ratio_failure(
    closure_context: SimpleNamespace,
    inner_open: bool,
) -> None:
    """Gamma1, alpha1 and A_down2 close at f1=1/6,r=1/2, not r=1.

    An intentionally retained f1=1/2 instead has inner net[1,-1]. The
    diagnostic must not call that an exclusively post-solve failure. Flows
    and sources are verified independently, with roundoff-only tolerances.
    """
    context = closure_context
    if inner_open:
        context.fractions = np.array([0.5, 0.5])
    evaluator = context.evaluator
    result = evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    saved_levels = result.absolute_level_population.copy(deep=True)
    saved_ions = result.ion_population.copy(deep=True)
    report = diagnose_level_closure(evaluator, 17)
    assert report["context"].loc[17, "branch"] == (
        "inner_level_solve_open" if inner_open else "post_level_state_changed"
    )
    assert bool(report["context"].loc[17, "terminal_replay_matches_reported"])
    np.testing.assert_allclose(
        report["levels"]["inner"]["net_flow"],
        [1.0, -1.0] if inner_open else [0.0, 0.0],
        atol=2e-16,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        report["levels"]["terminal"]["net_flow"]
        - report["levels"]["represented_fractions"]["net_flow"],
        [0.5, 0.5],
        atol=2e-16,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        report["stages"].loc["inner", "used_ion_ratio"], 0.5, rtol=1e-15
    )
    np.testing.assert_equal(
        report["stages"].loc["terminal", "used_ion_ratio"], 1.0
    )
    pdt.assert_frame_equal(
        result.absolute_level_population, saved_levels, check_exact=True
    )
    pdt.assert_frame_equal(result.ion_population, saved_ions, check_exact=True)
    assert evaluator.last_level_snapshot.evaluation is result


def test_diagnostic_uses_final_density_for_recombination(
    closure_context: SimpleNamespace,
) -> None:
    """Doubling ne doubles alpha; the eliminated ratio halves to preserve flux."""
    context = closure_context
    context.final_density = 2.0
    context.evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    report = diagnose_level_closure(context.evaluator, 17)
    np.testing.assert_allclose(
        report["rates"]["inner"]["recombination"], 1.0, rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        report["rates"]["final_density"]["recombination"],
        2.0,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        report["stages"].loc["final_density", "eliminated_ion_ratio"],
        0.25,
        rtol=1e-15,
    )
    np.testing.assert_allclose(
        report["levels"]["terminal"]["net_flow"], 1.5, rtol=1e-15
    )


def test_diagnostic_isolates_opacity_geometry(
    closure_context: SimpleNamespace,
) -> None:
    """Inner tau2 and final tau1.5 give independently known escape factors.

    f is prescribed for this sensitivity test, not claimed to solve the newly
    trapped atom. Geometry alone changes beta before the final ratio override.
    """
    context = closure_context
    context.bound = replace(
        context.bound,
        radiative_rate_coefficient=np.array([2.0]),
        fixed_rate=np.zeros(1),
        beta_line_idx=np.array([0]),
    )
    context.evaluator.sobolev_inputs = (
        replace(
            context.evaluator.sobolev_inputs[0], tau_coefficient=np.ones(1)
        ),
    )
    context.evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    report = diagnose_level_closure(context.evaluator, 17)
    np.testing.assert_allclose(
        report["beta"]["inner"], -np.expm1(-2.0) / 2.0, rtol=1e-14
    )
    np.testing.assert_allclose(
        report["beta"]["final_geometry"], -np.expm1(-1.5) / 1.5, rtol=1e-14
    )
    expected_bound_inflow = 2.0 * (-np.expm1(-1.5) / 1.5) / 6.0
    np.testing.assert_allclose(
        report["levels"]["final_geometry"]["net_flow"].iloc[0],
        expected_bound_inflow - 5 / 6 + 0.5,
        rtol=1e-13,
    )


def test_diagnostic_identifies_represented_fraction_change(
    closure_context: SimpleNamespace,
) -> None:
    """Exact f=2/3,1/3 becomes q/q; only reconstructed fractions open closure."""
    context = closure_context
    quantum = np.nextafter(0.0, 1.0)
    context.fractions = np.array([2 / 3, 1 / 3])
    context.ions.iloc[:, 0] = [2 * quantum, 0.0]
    context.bound = BoundBoundMatrixRates(
        2,
        np.array([0, 1]),
        np.array([1, 0]),
        np.zeros(2),
        np.zeros(2),
        np.array([1.0, 2.0]),
        np.array([-1, -1]),
    )
    context.continuum = replace(context.continuum, photoionization=np.zeros(2))
    context.evaluator.shell_number_densities = (
        ShellNumberDensity(2 * quantum, np.full(2, quantum), np.array([0, 1])),
    )
    context.evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    report = diagnose_level_closure(context.evaluator, 17)
    np.testing.assert_allclose(
        report["stages"].loc[
            ["inner", "final_density", "final_geometry"],
            "max_relative_residual",
        ],
        0.0,
        rtol=0.0,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        report["levels"]["represented_fractions"]["fraction"],
        0.5,
        atol=0.0,
        rtol=0.0,
    )
    allowance = 64 * np.finfo(float).eps / 1e-8 * 3.0
    np.testing.assert_allclose(
        report["stages"].loc["terminal", "max_relative_residual"],
        0.5 / (1.5 + allowance),
        rtol=1e-13,
    )


def test_diagnostic_default_never_calls_global_or_level_solvers(
    closure_context: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Read-only residual replay is distinct from expensive solving (orchestration)."""
    context = closure_context
    context.evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )

    def forbid_solve(*args: object, **kwargs: object) -> None:
        pytest.fail("Diagnostic unexpectedly invoked a solver")

    for name in [
        "evaluate",
        "_calculate_level_solution",
        "_solve_charge_population",
        "calculate_continuum_coefficients",
    ]:
        monkeypatch.setattr(context.evaluator, name, forbid_solve)
    report = diagnose_level_closure(context.evaluator, 17)
    assert "resolved_inner" not in report["stages"].index


def test_snapshot_replaces_only_completed_calls(
    closure_context: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep one completed call, not a history; a later exception preserves it."""
    context = closure_context
    evaluator = context.evaluator
    seed = pd.DataFrame(context.fractions, columns=[17])
    first = evaluator.evaluate(np.array([1.0]), np.array([1e4]), seed)
    second = evaluator.evaluate(np.array([1.0]), np.array([1e4]), seed)
    assert evaluator.last_level_snapshot.evaluation is second
    assert evaluator.last_level_snapshot.evaluation is not first
    assert evaluator.last_level_snapshot.evaluation_count == 2
    with pytest.raises(ValueError, match="positive"):
        evaluator.evaluate(np.array([-1.0]), np.array([1e4]), seed)
    assert evaluator.last_level_snapshot.evaluation is second
    assert (
        diagnose_level_closure(evaluator, 17)["context"].loc[
            17, "evaluation_count"
        ]
        == 2
    )


def test_diagnostic_resolves_only_requested_shell_and_keeps_snapshot(
    closure_context: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Retry only the requested position, without modifying the saved state.

    Duplicate the synthetic snapshot into labels91/17. A prescribed retry fixes
    the open inner distribution; its1/6 excited fraction closes the analytic
    atom. This is orchestration, not optimizer-performance validation.
    """
    context = closure_context
    context.fractions = np.array([0.5, 0.5])
    evaluator = context.evaluator
    result = evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    replacements = {}
    for name in [
        "normalized_population",
        "absolute_level_population",
        "ion_population",
        "relative_level_residual",
        "trial_beta_sobolev",
        "beta_sobolev",
    ]:
        frame = getattr(result, name)
        duplicated = pd.concat([frame, frame], axis=1)
        duplicated.columns = [91, 17]
        replacements[name] = duplicated
    replacements["trial_electron_density"] = pd.Series(
        [1.0, 1.0], index=[91, 17]
    )
    replacements["charge_solved_electron_density"] = pd.Series(
        [1.0, 1.0], index=[91, 17]
    )
    snapshot = evaluator.last_level_snapshot
    evaluator.last_level_snapshot = replace(
        snapshot,
        evaluation=replace(result, **replacements),
        electron_temperature=np.array([1e4, 1e4]),
        shell_number_densities=snapshot.shell_number_densities * 2,
        continuum_rates=snapshot.continuum_rates * 2,
        bound_bound_rates=snapshot.bound_bound_rates * 2,
    )
    evaluator.sobolev_inputs *= 2
    selected_positions = []

    def retry(*args: object) -> tuple:
        selected_positions.append(args[0])
        return np.array([5 / 6, 1 / 6]), 0.5, np.ones(1), np.zeros(2)

    monkeypatch.setattr(evaluator, "_calculate_level_solution", retry)
    retained = evaluator.last_level_snapshot
    report = diagnose_level_closure(evaluator, 17, resolve=True)
    assert selected_positions == [1]
    assert evaluator.last_level_snapshot is retained
    np.testing.assert_allclose(
        report["stages"].loc["resolved_inner", "max_relative_residual"],
        0.0,
        atol=2e-16,
    )
    assert (
        report["stages"].loc["resolved_terminal", "max_relative_residual"] > 0.1
    )
    assert report["context"].index.tolist() == [17]


def test_diagnostic_reports_stale_fixed_input_replay(
    closure_context: SimpleNamespace,
) -> None:
    """Changing recombination policy must invalidate diagnostic attribution."""
    context = closure_context
    context.evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    context.evaluator.case_b_species = frozenset([(1, 0)])
    report = diagnose_level_closure(context.evaluator, 17)
    assert report["context"].loc[17, "branch"] == "replay_context_mismatch"
    assert not bool(
        report["context"].loc[17, "terminal_replay_matches_reported"]
    )


def test_diagnostic_requires_a_completed_snapshot(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """A fresh evaluator cannot diagnose a historical run from logs alone."""
    with pytest.raises(ValueError, match="No completed"):
        diagnose_level_closure(toy_evaluator, 17)


def test_diagnostic_real_one_shell_retry_closes_analytic_inner_problem(
    closure_context: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The real level optimizer recovers f1=1/6, but cannot fix prescribed r=1.

    This verifies the optional diagnostic retry, not convergence of the global
    plasma. The independent two-level equation is -3*f1+1/2=0. Tolerance1e-10
    follows the existing stricter inner target and binary64 solve arithmetic.
    """
    context = closure_context
    context.fractions = np.array([0.5, 0.5])
    evaluator = context.evaluator
    evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    monkeypatch.setattr(
        evaluator,
        "_calculate_level_solution",
        PlasmaEquilibriumEvaluator._calculate_level_solution.__get__(evaluator),
    )
    report = diagnose_level_closure(evaluator, 17, resolve=True)
    np.testing.assert_allclose(
        report["levels"]["resolved_inner"]["fraction"],
        [5 / 6, 1 / 6],
        rtol=1e-10,
        atol=0.0,
    )
    assert (
        report["stages"].loc["resolved_inner", "max_relative_residual"]
        < evaluator.level_residual_tolerance
    )
    assert (
        report["stages"].loc["resolved_terminal", "max_relative_residual"] > 0.1
    )


def test_diagnostic_confirms_consistent_closed_state(
    closure_context: SimpleNamespace,
) -> None:
    """With authoritative ratio1/2, every staged physical equation closes."""
    context = closure_context
    context.ions.iloc[:, 0] = [8 / 3, 4 / 3]
    context.evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    report = diagnose_level_closure(context.evaluator, 17)
    assert report["context"].loc[17, "branch"] == "closed"
    assert bool(report["context"].loc[17, "terminal_replay_matches_reported"])
    np.testing.assert_allclose(
        report["stages"]["max_relative_residual"], 0.0, rtol=0.0, atol=2e-16
    )


def test_zero_selected_density_does_not_claim_authoritative_closure(
    closure_context: SimpleNamespace,
) -> None:
    """Expose legacy eliminated-ratio fallback when no finite actual ratio exists.

    The prescribed zero-neutral/positive-upper boundary is deliberately not an
    ionization solution. Do not mistake its zero legacy residual for verified
    authoritative balance; this test only verifies the diagnostic's warning.
    """
    context = closure_context
    context.ions.iloc[:, 0] = [0.0, 4.0]
    context.evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    report = diagnose_level_closure(context.evaluator, 17)
    assert report["context"].loc[17, "branch"] == "zero_selected_density"
    assert not bool(report["context"].loc[17, "authoritative_ratio_available"])
