"""Verify dimensional acceptance of selected-ion allocation underflow."""

from dataclasses import replace
from decimal import Decimal
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u

from tardis.plasma.equilibrium.evaluator import PlasmaEquilibriumEvaluator
from tardis.plasma.equilibrium.inputs import (
    BoundBoundMatrixRates,
    ContinuumRateCoefficients,
    ShellNumberDensity,
)
from tardis.plasma.equilibrium.solver import PlasmaEquilibriumSolver
from tardis.plasma.equilibrium.tests.test_evaluator import (
    toy_evaluator as toy_evaluator,  # noqa: PLC0414
)
from tardis.plasma.equilibrium.tests.test_level_failure_policy import (
    terminal_policy_evaluator as terminal_policy_evaluator,  # noqa: PLC0414
)


@pytest.fixture
def subnormal_evaluator(
    toy_evaluator: PlasmaEquilibriumEvaluator, monkeypatch: pytest.MonkeyPatch
) -> SimpleNamespace:
    """Provide an exact 1:2 atom and a prescribed authoritative ion boundary.

    The actual allocator/evaluator operates on two-q neutral density. The ion
    boundary supplies populations but no charge solution; charge validation is
    separately covered by the terminal solver's existing analytic fixtures.
    """
    evaluator = toy_evaluator
    quantum = np.nextafter(0.0, 1.0)
    fractions = np.array([2 / 3, 1 / 3])
    evaluator.case_b_species = frozenset()
    evaluator.time_explosion = 1e6 * u.s
    # This opt-in attribute is intentionally assigned after construction so
    # pre-change tests fail on scientific acceptance, not a missing argument.
    evaluator.absolute_level_population_tolerance = 1e-20
    evaluator.sobolev_inputs = (
        replace(evaluator.sobolev_inputs[0], tau_coefficient=np.zeros(1)),
    )
    evaluator.shell_number_densities = (
        ShellNumberDensity(
            2 * quantum, 2 * quantum * fractions, np.array([0, 1])
        ),
    )
    evaluator.elemental_number_density.iloc[0, 0] = 2 * quantum
    context = SimpleNamespace(
        evaluator=evaluator,
        seed=pd.DataFrame(fractions, index=evaluator.level_population_index),
        ions=pd.DataFrame(
            [[2 * quantum], [0.0]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0), (1, 1)], names=["atomic_number", "ion_number"]
            ),
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
        continuum=ContinuumRateCoefficients(
            np.zeros(2), np.zeros(2), np.ones(2), np.zeros(2), np.zeros(2)
        ),
    )
    original_coefficients = evaluator.calculate_continuum_coefficients

    def fixed_continuum(temperature: npt.NDArray[np.float64]) -> tuple:
        coefficients = original_coefficients(temperature)
        return ((context.continuum,), *coefficients[1:])

    monkeypatch.setattr(
        evaluator, "calculate_continuum_coefficients", fixed_continuum
    )
    monkeypatch.setattr(
        evaluator.rate_matrix_solver,
        "prepare_species_rates",
        lambda *args: (context.bound,),
    )
    monkeypatch.setattr(
        evaluator,
        "_solve_charge_population",
        lambda *args: (context.ions.copy(deep=True), None),
    )
    return context


def test_subnormal_absolute_budget_retains_relative_defect(
    subnormal_evaluator: SimpleNamespace,
) -> None:
    """Parent2q cannot have both exact conservation and 1:2 level balance.

    a+b=2 and a=2b have no integer-quanta solution. The allocator gives q/q;
    net/gross is +/-1/3 before the existing rate-roundoff allowance. Accept
    this negligible absolute turnover, without changing densities or residuals.
    """
    context = subnormal_evaluator
    evaluator = context.evaluator
    result = evaluator.evaluate(np.array([1.0]), np.array([1e4]), context.seed)
    quantum = np.nextafter(0.0, 1.0)
    np.testing.assert_array_equal(result.absolute_level_population, quantum)
    denominator = 3.0 + (64 * np.finfo(float).eps / 1e-8) * 6.0
    np.testing.assert_allclose(
        result.relative_level_residual.to_numpy().ravel(),
        [1 / denominator, -1 / denominator],
        rtol=1e-13,
        atol=0.0,
    )
    budget = getattr(result, "absolute_level_budget", None)
    assert budget is not None and bool(budget.loc[0, "accepted"])
    pdt.assert_frame_equal(
        result.ion_population, context.ions, check_exact=True
    )


@pytest.mark.parametrize("enabled", [False, True])
@pytest.mark.parametrize("age", [None, 0.0, -1.0, np.inf, 1e6])
def test_absolute_budget_requires_explicit_tolerance_and_usable_age(
    subnormal_evaluator: SimpleNamespace, enabled: bool, age: float | None
) -> None:
    """Generic/default and missing-time contexts never gain an implicit floor."""
    context = subnormal_evaluator
    context.evaluator.absolute_level_population_tolerance = (
        1e-20 if enabled else 0.0
    )
    context.evaluator.time_explosion = None if age is None else age * u.s
    result = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    budget = getattr(result, "absolute_level_budget", None)
    accepted = budget is not None and bool(budget.loc[0, "accepted"])
    assert accepted is (enabled and age == 1e6)


@pytest.mark.parametrize("source", [False, True])
def test_absolute_budget_rejects_significant_turnover(
    subnormal_evaluator: SimpleNamespace, source: bool
) -> None:
    """Subnormal population is not permission to ignore relevant absolute flow.

    Without a source, q*0.25 rounds to zero but its integrated turnover over
    1e308 s is ~1e-16 cm^-3, larger than1e-20. With a source, recombination
    from the prescribed normal upper-ion density dominates even at zero lower
    ionization. Decimal-free logarithmic oracles avoid premature products.
    """
    context = subnormal_evaluator
    context.evaluator.time_explosion = (1e200 if source else 1e308) * u.s
    context.bound = replace(context.bound, fixed_rate=np.array([0.25, 0.5]))
    if source:
        context.ions.iloc[1, 0] = 1e-200
    result = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    budget = getattr(result, "absolute_level_budget", None)
    assert budget is not None
    assert not bool(budget.loc[0, "accepted"])
    assert budget.loc[0, "log_integrated_turnover"] > np.log(1e-20)


def test_absolute_budget_rejects_non_cancelled_opacity(
    subnormal_evaluator: SimpleNamespace,
) -> None:
    """A subnormal atom with finite coefficient1e308 can have order-unity tau.

    The stimulated subtraction need not be large for its separate absolute
    contributions to be significant. No numerical density-only skip is safe.
    """
    context = subnormal_evaluator
    parent = 0.75 * np.finfo(float).tiny
    context.seed.iloc[:, 0] = [1 / 3, 2 / 3]
    context.bound = replace(context.bound, fixed_rate=np.array([2.0, 1.0]))
    context.ions.iloc[0, 0] = parent
    context.evaluator.shell_number_densities = (
        ShellNumberDensity(
            parent, parent * np.array([1 / 3, 2 / 3]), np.array([0, 1])
        ),
    )
    context.evaluator.elemental_number_density.iloc[0, 0] = parent
    context.evaluator.sobolev_inputs = (
        replace(
            context.evaluator.sobolev_inputs[0],
            tau_coefficient=np.array([1e308]),
        ),
    )
    result = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    budget = getattr(result, "absolute_level_budget", None)
    assert budget is not None and not bool(budget.loc[0, "accepted"])
    assert budget.loc[0, "log_optical_depth_bound"] > np.log(1e-8)
    np.testing.assert_equal(result.tau_sobolev.iloc[0, 0], 0.0)


@pytest.mark.parametrize("parent", [np.finfo(float).tiny, 1.0])
def test_absolute_budget_never_excuses_normal_parent(
    subnormal_evaluator: SimpleNamespace, parent: float
) -> None:
    """Eligibility concerns the entire parent, not just individual weak levels."""
    context = subnormal_evaluator
    context.ions.iloc[0, 0] = parent
    result = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    budget = getattr(result, "absolute_level_budget", None)
    assert budget is not None and not bool(budget.loc[0, "accepted"])


@pytest.mark.parametrize("fraction_of_budget", [0.99, 1.01])
def test_absolute_budget_uses_dimensional_turnover_boundary(
    subnormal_evaluator: SimpleNamespace, fraction_of_budget: float
) -> None:
    """G=6q s^-1 integrated over1e308s is compared before products vanish.

    Decimal computes the independent physical product from exact binary64 q;
    a1% separation is far larger than log-domain arithmetic uncertainty.
    """
    context = subnormal_evaluator
    context.evaluator.time_explosion = 1e308 * u.s
    exact_turnover = (
        Decimal.from_float(np.nextafter(0.0, 1.0))
        * 6
        * Decimal.from_float(1e308)
    )
    context.evaluator.absolute_level_population_tolerance = (
        float(exact_turnover) / fraction_of_budget
    )
    result = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    budget = result.absolute_level_budget
    np.testing.assert_allclose(
        budget.loc[0, "log_integrated_turnover"],
        np.log(float(exact_turnover)),
        rtol=0.0,
        atol=5e-13,
    )
    assert bool(budget.loc[0, "accepted"]) is (fraction_of_budget < 1.0)


def test_absolute_budget_preserves_primitive_recombination_products(
    subnormal_evaluator: SimpleNamespace,
) -> None:
    """A coefficientq times ne=.25 underflows, but its source is significant.

    q*.25*1e308*1e6 is ~1.235e-10 cm^-3, exceeding1e-20. The direct
    component test isolates this budget from the pre-existing fractional
    residual's distinct continuum-ratio overflow limitations.
    """
    context = subnormal_evaluator
    snapshot = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    quantum = np.nextafter(0.0, 1.0)
    context.ions.iloc[1, 0] = 1e308
    continuum = replace(
        context.continuum, spontaneous_recombination=np.array([quantum, 0.0])
    )
    budget = context.evaluator._calculate_absolute_level_budget(
        pd.Series([0.25]),
        np.array([1e4]),
        snapshot.absolute_level_population,
        context.ions,
        SimpleNamespace(species_index=snapshot.normalized_population.index),
        SimpleNamespace(continuum_rate_coefficients=(continuum,)),
        snapshot.beta_sobolev,
    )
    assert not bool(budget.loc[0, "accepted"])
    oracle = float(
        Decimal.from_float(quantum)
        * Decimal("0.25")
        * Decimal("1e308")
        * Decimal("1e6")
    )
    np.testing.assert_allclose(
        budget.loc[0, "log_integrated_turnover"],
        np.log(oracle),
        rtol=0.0,
        atol=5e-13,
    )


@pytest.mark.parametrize("invalid", [-1.0, np.nan, np.inf])
def test_absolute_budget_does_not_excuse_invalid_primitive_rates(
    subnormal_evaluator: SimpleNamespace, invalid: float
) -> None:
    """A zero population cannot sanitize invalid source or transition inputs."""
    context = subnormal_evaluator
    snapshot = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    continuum = replace(
        context.continuum, spontaneous_recombination=np.array([invalid, 0.0])
    )
    budget = context.evaluator._calculate_absolute_level_budget(
        pd.Series([1.0]),
        np.array([1e4]),
        snapshot.absolute_level_population,
        context.ions,
        SimpleNamespace(species_index=snapshot.normalized_population.index),
        SimpleNamespace(continuum_rate_coefficients=(continuum,)),
        snapshot.beta_sobolev,
    )
    assert not bool(budget.loc[0, "accepted"])


def test_absolute_budget_recomputed_after_timescale_change(
    subnormal_evaluator: SimpleNamespace,
) -> None:
    """Acceptance belongs to the current shell context, not a cached mask."""
    context = subnormal_evaluator
    first = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    context.evaluator.time_explosion = 1e308 * u.s
    second = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    assert bool(first.absolute_level_budget.loc[0, "accepted"])
    assert not bool(second.absolute_level_budget.loc[0, "accepted"])


def test_absolute_budget_checks_represented_sum_and_shell_labels(
    subnormal_evaluator: SimpleNamespace,
) -> None:
    """Each subnormal level is insufficient if their sum is a normal density.

    Nonconsecutive shell labels keep an eligible first shell separate from
    a second one with deliberately inconsistent normal represented sum. The
    terminal population guards independently reject that inconsistency.
    """
    context = subnormal_evaluator
    snapshot = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    levels = pd.concat(
        [
            snapshot.absolute_level_population,
            snapshot.absolute_level_population,
        ],
        axis=1,
    )
    levels.columns = [91, 54]
    levels.loc[:, 54] = 0.5 * np.finfo(float).tiny
    ions = pd.concat([context.ions, context.ions], axis=1)
    ions.columns = levels.columns
    context.evaluator.sobolev_inputs *= 2
    budget = context.evaluator._calculate_absolute_level_budget(
        pd.Series([1.0, 1.0], index=levels.columns),
        np.array([1e4, 1e4]),
        levels,
        ions,
        SimpleNamespace(species_index=snapshot.normalized_population.index),
        SimpleNamespace(continuum_rate_coefficients=(context.continuum,) * 2),
        pd.DataFrame([[1.0, 1.0]], columns=levels.columns),
    )
    pd.testing.assert_index_equal(budget.index, levels.columns)
    np.testing.assert_array_equal(budget["accepted"], [True, False])


@pytest.mark.parametrize(
    "corruption", ["missing_boolean", "wrong_shells", "nonfinite_residual"]
)
def test_absolute_budget_cannot_hide_invalid_acceptance_data(
    corruption: str,
) -> None:
    """NaN flags, misaligned labels and nonfinite residuals cannot become true."""
    snapshot = SimpleNamespace(
        relative_level_residual=pd.DataFrame([[0.3, 0.3]], columns=[91, 54]),
        absolute_level_budget=pd.DataFrame(
            {"accepted": [True, False]}, index=[91, 54]
        ),
    )
    if corruption == "missing_boolean":
        snapshot.absolute_level_budget = pd.DataFrame(
            {"accepted": [np.nan, False]}, index=[91, 54]
        )
    elif corruption == "wrong_shells":
        snapshot.absolute_level_budget.index = [54, 91]
    else:
        snapshot.relative_level_residual.iloc[0, 0] = np.nan
    if corruption == "nonfinite_residual":
        np.testing.assert_array_equal(
            PlasmaEquilibriumSolver._calculate_level_acceptance(snapshot, 1e-8),
            [[False, False]],
        )
    else:
        with pytest.raises(RuntimeError, match="Absolute level budget"):
            PlasmaEquilibriumSolver._calculate_level_acceptance(snapshot, 1e-8)


def test_terminal_absolute_acceptance_skips_refinement(
    subnormal_evaluator: SimpleNamespace, terminal_policy_evaluator: object
) -> None:
    """Actual absolute budget skips retries without overwriting raw diagnostics.

    The analytic-charge boundary isolates solver orchestration; physical
    derivation of its supplied acceptance is checked separately above.
    """
    context = subnormal_evaluator
    snapshot = context.evaluator.evaluate(
        np.array([1.0]), np.array([1e4]), context.seed
    )
    actual_budget = getattr(snapshot, "absolute_level_budget", None)
    evaluator = terminal_policy_evaluator
    original_evaluate = evaluator.evaluate

    def evaluate_with_budget(*args: object) -> SimpleNamespace:
        result = original_evaluate(*args)
        if actual_budget is not None:
            result.absolute_level_budget = actual_budget.set_axis([17])
        return result

    evaluator.evaluate = evaluate_with_budget
    result = PlasmaEquilibriumSolver(fixed_temperature_ratio=0.9).solve(
        evaluator,
        [4e8],
        [9000.0],
        pd.DataFrame([[0.75], [0.25]], columns=[17]),
        [2e9],
        [1e4],
    )
    assert evaluator.refinements == 0
    assert result.level_converged
    np.testing.assert_equal(result.maximum_relative_level_residual, 0.1921)
