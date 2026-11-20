"""Verify seed reuse without caching acceptance of changed kinetic equations."""

from collections.abc import Callable
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from scipy.optimize import OptimizeResult, brentq
from scipy.optimize import root as scipy_root

from tardis.plasma.equilibrium.diagnostics import diagnose_level_closure
from tardis.plasma.equilibrium.evaluator import PlasmaEquilibriumEvaluator
from tardis.plasma.equilibrium.inputs import (
    BoundBoundMatrixRates,
    ContinuumRateCoefficients,
)
from tardis.plasma.equilibrium.tests.test_closure_diagnostic import (
    closure_context as closure_context,  # noqa: PLC0414
)
from tardis.plasma.equilibrium.tests.test_evaluator import (
    toy_evaluator as toy_evaluator,  # noqa: PLC0414
)


@pytest.fixture
def warm_atom(toy_evaluator: PlasmaEquilibriumEvaluator) -> SimpleNamespace:
    """Provide a two-level optically thin atom with exact up/down balance1/2."""
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()
    return SimpleNamespace(
        evaluator=evaluator,
        continuum=ContinuumRateCoefficients(
            np.zeros(2), np.zeros(2), np.ones(2), np.zeros(2), np.zeros(2)
        ),
        bound=BoundBoundMatrixRates(
            2,
            np.array([0, 1]),
            np.array([1, 0]),
            np.zeros(2),
            np.zeros(2),
            np.array([1.0, 2.0]),
            np.full(2, -1),
        ),
        geometry=evaluator.shell_number_densities[0],
        sobolev=replace(
            evaluator.sobolev_inputs[0], tau_coefficient=np.zeros(1)
        ),
    )


def test_repeated_closed_atom_reuses_seed_without_linear_initialization(
    warm_atom: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A repeated atom retains exact f1=1/3 rather than resolve a poor seed.

    The physical oracle is upward f0 = downward2*f1 and sum(f)=1.
    Population rtol1e-8 follows the current kinetic tolerance. Linear-solve
    counts separately verify avoiding redundant initialization, not wall time.
    """
    context = warm_atom
    original = np.linalg.solve
    calls = []

    def count_linear(
        *args: object, **kwargs: object
    ) -> npt.NDArray[np.float64]:
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(np.linalg, "solve", count_linear)
    seed = np.array([0.9, 0.1])
    arguments = (
        0,
        1.0,
        1e4,
        seed,
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )
    first = context.evaluator._calculate_level_solution(*arguments)[0]
    assert len(calls) == 1
    calls.clear()
    second = context.evaluator._calculate_level_solution(*arguments)[0]
    assert_allclose(first, [2 / 3, 1 / 3], rtol=1e-8, atol=0.0)
    assert_allclose(second, [2 / 3, 1 / 3], rtol=1e-8, atol=0.0)
    assert_allclose(seed, [0.9, 0.1], rtol=0.0, atol=0.0)
    assert not calls, (
        "Repeated closed atom unnecessarily rebuilt the frozen-beta initializer"
    )


@pytest.mark.parametrize("change", ["rates", "density", "opacity"])
def test_warm_seed_rechecks_changed_equations(
    warm_atom: SimpleNamespace,
    change: str,
) -> None:
    """Changed rates/density/opacity require a new equilibrium, not a cache hit.

    Finite continuum has exact eliminated scalar balance; nonlinear trapping
    is independently solved with brentq and analytic beta. rtol1e-8 matches
    the existing physical tolerance, not an added acceptance floor.
    """
    context = warm_atom
    evaluator = context.evaluator
    seed = np.array([0.9, 0.1])
    evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        seed,
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )
    density = 1.0
    expected = np.array([0.4, 0.6])
    if change == "rates":
        context.bound = replace(context.bound, fixed_rate=np.array([3.0, 2.0]))
    elif change == "density":
        density = 3.0
        context.continuum = ContinuumRateCoefficients(
            np.array([2.0, 4.0]),
            np.array([0.5, 0.25]),
            np.array([1.0, 2.0]),
            np.zeros(2),
            np.zeros(2),
        )
        gamma = np.array([2.0, 4.0]) + density * np.array([0.5, 0.25])
        excited = (1.0 + 2 / 3 * gamma[0]) / (
            3.0 + gamma[1] - 2 / 3 * (gamma[1] - gamma[0])
        )
        expected = np.array([1 - excited, excited])
    else:
        context.bound = replace(
            context.bound,
            fixed_rate=np.array([1.0, 0.0]),
            radiative_rate_coefficient=np.array([0.0, 2.0]),
            beta_line_idx=np.array([-1, 0]),
        )
        context.sobolev = replace(
            context.sobolev, tau_coefficient=np.array([1e-9])
        )

        def scalar_balance(excited: float) -> float:
            tau = 10 * (1 - 1.5 * excited)
            beta = -np.expm1(-tau) / tau if tau != 0 else 1.0
            return 1 - excited - 2 * beta * excited

        excited = brentq(scalar_balance, 0.0, 2 / 3, xtol=1e-14)
        expected = np.array([1 - excited, excited])
    warm = evaluator._calculate_level_solution(
        0,
        density,
        1e4,
        seed,
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )[0]
    assert_allclose(warm, expected, rtol=1e-8, atol=0.0)
    evaluator._level_solution_cache.clear()
    cold = evaluator._calculate_level_solution(
        0,
        density,
        1e4,
        seed,
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )[0]
    assert_allclose(cold, expected, rtol=1e-8, atol=0.0)
    assert_allclose(warm, cold, rtol=1e-8, atol=0.0)


def test_better_caller_seed_survives_stale_cache(
    warm_atom: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exact caller f=[.4,.6] outranks a formerly correct1/3 excited seed."""
    context = warm_atom
    context.evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        np.array([0.9, 0.1]),
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )
    bound = replace(context.bound, fixed_rate=np.array([3.0, 2.0]))

    def forbid_initializer(*args: object, **kwargs: object) -> None:
        pytest.fail("A closed caller seed was displaced by a worse cached seed")

    monkeypatch.setattr(np.linalg, "solve", forbid_initializer)
    actual = context.evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        np.array([0.4, 0.6]),
        context.continuum,
        context.geometry,
        context.sobolev,
        bound,
    )[0]
    assert_allclose(actual, [0.4, 0.6], rtol=1e-13, atol=0.0)


def test_warm_seed_enters_nonlinear_fallback_when_initializer_unavailable(
    warm_atom: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A better cached trial actually seeds HYBR after a frozen-linear failure.

    Former f1=1/3 is closer than caller f1=.1 to the new exact equilibrium.6.
    Root-input counting verifies warm-seed orchestration; the analytic final
    distribution independently verifies the unchanged physical equations.
    """
    context = warm_atom
    evaluator = context.evaluator
    seed = np.array([0.9, 0.1])
    evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        seed,
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )
    bound = replace(context.bound, fixed_rate=np.array([3.0, 2.0]))
    root_seeds = []

    def unavailable_initializer(*args: object, **kwargs: object) -> None:
        raise np.linalg.LinAlgError("controlled optional initializer failure")

    def record_root(
        function: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
        initial: npt.NDArray[np.float64],
        **kwargs: object,
    ) -> OptimizeResult:
        root_seeds.append(initial.copy())
        return scipy_root(function, initial, **kwargs)

    monkeypatch.setattr(np.linalg, "solve", unavailable_initializer)
    monkeypatch.setattr("tardis.plasma.equilibrium.evaluator.root", record_root)
    result = evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        seed,
        context.continuum,
        context.geometry,
        context.sobolev,
        bound,
    )[0]
    assert_allclose(root_seeds[0], [2 / 3, 1 / 3], rtol=1e-13, atol=0.0)
    assert_allclose(result, [0.4, 0.6], rtol=1e-8, atol=0.0)


def test_frozen_initializer_uses_opacity_of_best_current_seed(
    warm_atom: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ranking a worse cache last must not freeze its stale beta.

    At caller f=[.45,.55], tau=10*(.45-.55/2)=1.75. The independently
    row-scaled second matrix row is [1,-2*beta(1.75)], not the cache's beta.
    This catches a last-evaluation versus best-seed bookkeeping mistake.
    """
    context = warm_atom
    evaluator = context.evaluator
    evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        np.array([0.9, 0.1]),
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )
    bound = replace(
        context.bound,
        fixed_rate=np.array([1.0, 0.0]),
        radiative_rate_coefficient=np.array([0.0, 2.0]),
        beta_line_idx=np.array([-1, 0]),
    )
    sobolev = replace(context.sobolev, tau_coefficient=np.array([1e-9]))
    matrices = []
    original = np.linalg.solve

    def record_linear(
        matrix: npt.NDArray[np.float64], rhs: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        matrices.append(matrix.copy())
        return original(matrix, rhs)

    monkeypatch.setattr(np.linalg, "solve", record_linear)
    evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        np.array([0.45, 0.55]),
        context.continuum,
        context.geometry,
        sobolev,
        bound,
    )
    beta = -np.expm1(-1.75) / 1.75
    assert_allclose(matrices[0][1], [1.0, -2 * beta], rtol=1e-13, atol=0.0)


def test_cache_uses_shell_labels_and_prunes_removed_shells(
    warm_atom: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reordered positions retain their own exact atom, not another shell's.

    Shell91 has up/down1/2 and shell17 up/down3/2. Distinct exact fractions
    and a forbidden initializer detect accidental positional reuse.
    """
    context = warm_atom
    evaluator = context.evaluator
    evaluator.elemental_number_density = pd.DataFrame(
        [[1e10, 1e10]], index=[1], columns=[91, 17]
    )
    second_bound = replace(context.bound, fixed_rate=np.array([3.0, 2.0]))
    for position, bound in enumerate([context.bound, second_bound]):
        evaluator._calculate_level_solution(
            position,
            1.0,
            1e4,
            np.array([0.9, 0.1]),
            context.continuum,
            context.geometry,
            context.sobolev,
            bound,
        )
    evaluator.elemental_number_density = evaluator.elemental_number_density[
        [17, 91]
    ]

    def forbid_initializer(*args: object, **kwargs: object) -> None:
        pytest.fail("Reordering lost the successful seed for this shell label")

    monkeypatch.setattr(np.linalg, "solve", forbid_initializer)
    reordered = evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        np.array([0.9, 0.1]),
        context.continuum,
        context.geometry,
        context.sobolev,
        second_bound,
    )[0]
    assert_allclose(reordered, [0.4, 0.6], rtol=1e-8, atol=0.0)
    evaluator.elemental_number_density = (
        evaluator.elemental_number_density.reindex(
            columns=[17, 5], fill_value=1e10
        )
    )
    evaluator._calculate_level_solution(
        0,
        1.0,
        1e4,
        np.array([0.9, 0.1]),
        context.continuum,
        context.geometry,
        context.sobolev,
        second_bound,
    )
    assert set(evaluator._level_solution_cache) == {17}


@pytest.mark.parametrize("change", ["species", "level_order"])
def test_cache_invalidates_structural_context(
    warm_atom: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
    change: str,
) -> None:
    """Identity changes must not reuse unlabeled old fractions (orchestration)."""
    context = warm_atom
    evaluator = context.evaluator
    arguments = (
        0,
        1.0,
        1e4,
        np.array([0.9, 0.1]),
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )
    evaluator._calculate_level_solution(*arguments)
    if change == "species":
        # The rate helper accepts supplied arrays independently of atomic Z;
        # this isolates cache context, not a valid all-data plasma rebuild.
        evaluator.equilibrium_species = (2, 0)
    else:
        evaluator.level_population_index = evaluator.level_population_index[
            ::-1
        ]
    original = np.linalg.solve
    calls = []

    def count_initializer(
        *args: object, **kwargs: object
    ) -> npt.NDArray[np.float64]:
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(np.linalg, "solve", count_initializer)
    evaluator._calculate_level_solution(*arguments)
    assert len(calls) == 1
    assert len(evaluator._level_solution_cache) == 1


def test_cache_owns_fractions_independent_of_returned_array(
    warm_atom: SimpleNamespace,
) -> None:
    """Returned mutable populations cannot corrupt the next call's warm seed."""
    context = warm_atom
    seed = np.array([0.9, 0.1])
    arguments = (
        0,
        1.0,
        1e4,
        seed,
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )
    result = context.evaluator._calculate_level_solution(*arguments)[0]
    cached = context.evaluator._level_solution_cache[0]
    assert not np.shares_memory(cached, result)
    assert not np.shares_memory(cached, seed)
    result[:] = 0.0
    repeated = context.evaluator._calculate_level_solution(*arguments)[0]
    assert_allclose(repeated, [2 / 3, 1 / 3], rtol=1e-8, atol=0.0)


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("nonfinite", [False, True])
def test_open_or_nonfinite_candidate_does_not_poison_cache(
    warm_atom: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
    existing: bool,
    nonfinite: bool,
) -> None:
    """Only current all-row closure may store a warm seed, not solver success.

    Controlled optimizer and residual boundaries return a feasible but open
    candidate or no finite candidate. This tests cache admission/failure
    orchestration, not the scientific adequacy of those controlled rates.
    """
    context = warm_atom
    evaluator = context.evaluator
    arguments = (
        0,
        1.0,
        1e4,
        np.array([0.9, 0.1]),
        context.continuum,
        context.geometry,
        context.sobolev,
        context.bound,
    )
    if existing:
        evaluator._calculate_level_solution(*arguments)
    before = {
        key: value.copy()
        for key, value in evaluator._level_solution_cache.items()
    }

    def controlled_residual(
        *args: object, **kwargs: object
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float]:
        budget = kwargs.get("diagnostic_budget")
        values = np.full(2, np.nan) if nonfinite else np.array([0.1, 0.0])
        if budget is not None:
            budget.update(
                net_flow=values.copy(),
                gross_flow=np.ones(2),
                flow_scale=np.ones(2),
            )
        if not nonfinite and kwargs.get("replace_normalization_row", True):
            values[0] = 0.0
        return values, np.ones(1), 0.0

    monkeypatch.setattr(
        evaluator, "_calculate_level_state", controlled_residual
    )
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.root",
        lambda *args, **kwargs: SimpleNamespace(x=np.array([0.5, 0.5])),
    )
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.least_squares",
        lambda *args, **kwargs: SimpleNamespace(x=np.zeros(1)),
    )
    if nonfinite:
        with pytest.raises(ValueError, match="no physical candidate"):
            evaluator._calculate_level_solution(*arguments)
    else:
        result = evaluator._calculate_level_solution(*arguments)[0]
        assert_allclose(result.sum(), 1.0, rtol=0.0, atol=2e-16)
    assert set(evaluator._level_solution_cache) == set(before)
    for key, expected in before.items():
        assert_allclose(
            evaluator._level_solution_cache[key], expected, rtol=0.0, atol=0.0
        )


@pytest.mark.parametrize("raises", [False, True])
@pytest.mark.parametrize("present", [False, True])
def test_diagnostic_retry_restores_cache_and_context_even_on_failure(
    closure_context: SimpleNamespace,
    monkeypatch: pytest.MonkeyPatch,
    raises: bool,
    present: bool,
) -> None:
    """Read-only diagnostics cannot warm or invalidate subsequent plasma calls."""
    context = closure_context
    evaluator = context.evaluator
    evaluator.evaluate(
        np.array([1.0]),
        np.array([1e4]),
        pd.DataFrame(context.fractions, columns=[17]),
    )
    names = (
        "_level_solution_cache",
        "_level_solution_cache_species",
        "_level_solution_cache_index",
        "_level_solution_cache_columns",
    )
    if present:
        evaluator._level_solution_cache[17] = np.array([0.75, 0.25])
    else:
        for name in names:
            delattr(evaluator, name)
    saved = {
        name: getattr(evaluator, name)
        for name in names
        if hasattr(evaluator, name)
    }
    original_cache = saved.get("_level_solution_cache")
    original_value = None if original_cache is None else original_cache[17]

    def mutate_retry(
        *args: object,
    ) -> tuple[
        npt.NDArray[np.float64],
        float,
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        if present:
            original_cache[17][:] = [0.1, 0.9]
            original_cache[91] = np.ones(2)
        evaluator._level_solution_cache = {91: np.ones(2)}
        evaluator._level_solution_cache_species = (99, 0)
        evaluator._level_solution_cache_index = pd.Index([999])
        evaluator._level_solution_cache_columns = pd.Index([999])
        if raises:
            raise RuntimeError("controlled diagnostic retry failure")
        return context.fractions.copy(), 0.5, np.ones(1), np.zeros(2)

    monkeypatch.setattr(evaluator, "_calculate_level_solution", mutate_retry)
    if raises:
        with pytest.raises(RuntimeError, match="controlled diagnostic"):
            diagnose_level_closure(evaluator, 17, resolve=True)
    else:
        diagnose_level_closure(evaluator, 17, resolve=True)
    for name in names:
        if present:
            assert getattr(evaluator, name) is saved[name]
        else:
            assert not hasattr(evaluator, name)
    if present:
        assert set(original_cache) == {17}
        assert original_cache[17] is original_value
        assert_allclose(original_value, [0.75, 0.25], rtol=0.0, atol=0.0)
