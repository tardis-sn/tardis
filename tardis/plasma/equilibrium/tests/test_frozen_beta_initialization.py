"""Verify frozen-opacity seeding without replacing nonlinear level closure."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from scipy.optimize import brentq

from tardis.plasma.equilibrium.evaluator import PlasmaEquilibriumEvaluator
from tardis.plasma.equilibrium.inputs import (
    BoundBoundMatrixRates,
    ContinuumRateCoefficients,
    FixedRateState,
)
from tardis.plasma.equilibrium.tests.test_evaluator import (
    toy_evaluator as toy_evaluator,  # noqa: PLC0414 - shared pytest fixture
)


@pytest.mark.parametrize("rate_scale", [1e-12, 1.0, 1e12])
def test_frozen_beta_recovers_saturated_trace_level(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    rate_scale: float,
) -> None:
    """Exact up/down balance recovers a trace level from a saturated seed.

    Regime: two optically thin levels, up/down=1e-20/1 s^-1 and no populated
    continuum. Analytic f1=a/(1+a) and explicit net/gross fluxes independently
    verify the solution. Population rtol=1e-8 follows the configured physical
    tolerance; the same solution must survive uniform rate scaling.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()
    upward_rate = 1e-20 * rate_scale
    downward_rate = rate_scale
    bound_bound = BoundBoundMatrixRates(
        2,
        np.array([0, 1]),
        np.array([1, 0]),
        np.zeros(2),
        np.zeros(2),
        np.array([upward_rate, downward_rate]),
        np.full(2, -1),
    )
    continuum = ContinuumRateCoefficients(
        np.zeros(2),
        np.zeros(2),
        np.ones(2) * rate_scale,
        np.zeros(2),
        np.zeros(2),
    )
    seed = np.array([0.9999, 0.0001])
    fractions, _, _, _ = evaluator._calculate_level_solution(
        0,
        1e9,
        1e4,
        seed,
        continuum,
        evaluator.shell_number_densities[0],
        replace(evaluator.sobolev_inputs[0], tau_coefficient=np.zeros(1)),
        bound_bound,
    )
    inflow = upward_rate * fractions[0]
    outflow = downward_rate * fractions[1]
    gross = inflow + outflow
    allowance = 64 * np.finfo(float).eps / evaluator.level_residual_tolerance
    relative_residual = abs(inflow - outflow) / (gross * (1 + 2 * allowance))
    assert relative_residual <= evaluator.level_residual_tolerance
    assert_allclose(fractions[1], 1e-20 / (1 + 1e-20), rtol=1e-8, atol=0.0)
    assert_allclose(fractions.sum(), 1.0, rtol=0.0, atol=3e-16)
    assert np.all(fractions >= 0.0)
    assert_allclose(seed, [0.9999, 0.0001], rtol=0.0, atol=0.0)


@pytest.mark.parametrize("case_b", [False, True])
@pytest.mark.parametrize("electron_density", [2.0, 4.0])
def test_frozen_beta_seed_preserves_continuum_and_fixed_rates(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    case_b: bool,
    electron_density: float,
) -> None:
    """Solve independent bound/continuum balances with all density factors.

    Regime: two levels, bound up/down=2/6 s^-1 (including fixed excitation
    0.5), photoionization [2,4], fixed ionization 3 s^-1 and the displayed
    electron coefficients. Eliminating the third population analytically
    gives the excited fraction below. rtol=1e-12 allows linear-solve roundoff.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset({(1, 0)} if case_b else ())
    fixed_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
        ],
    )
    evaluator.fixed_rates = FixedRateState(
        ionization_rates=pd.DataFrame([[3.0]], index=fixed_index, columns=[0])
    )
    continuum = ContinuumRateCoefficients(
        np.array([2.0, 4.0]),
        np.array([0.5, 0.25]),
        np.array([1.0, 2.0]),
        np.array([0.1, 0.2]),
        np.array([0.2, 0.3]),
    )
    bound_bound = BoundBoundMatrixRates(
        2,
        np.array([0, 1]),
        np.array([1, 0]),
        np.zeros(2),
        np.array([1.5, 6.0]),
        np.array([0.5, 0.0]),
        np.full(2, -1),
    )
    gamma = (
        np.array([2.0 + 0.5 * electron_density, 4.0 + 0.25 * electron_density])
        + 3.0
    )
    alpha = electron_density * np.array(
        [1.1 + 0.2 * electron_density, 2.2 + 0.3 * electron_density]
    )
    if case_b:
        alpha[0] = 0.0
    expected_excited = (2.0 + alpha[1] * gamma[0] / alpha.sum()) / (
        8.0 + gamma[1] - alpha[1] * (gamma[1] - gamma[0]) / alpha.sum()
    )
    expected = np.array([1.0 - expected_excited, expected_excited])
    expected_ratio = gamma @ expected / alpha.sum()
    fractions, ratio, _, _ = evaluator._calculate_level_solution(
        0,
        electron_density,
        1e4,
        np.array([0.99, 0.01]),
        continuum,
        evaluator.shell_number_densities[0],
        replace(evaluator.sobolev_inputs[0], tau_coefficient=np.zeros(1)),
        bound_bound,
    )
    assert_allclose(fractions, expected, rtol=1e-12, atol=0.0)
    assert_allclose(ratio, expected_ratio, rtol=1e-12, atol=0.0)
    # The explicit continuum equation is also satisfied, not the seed's ratio.
    assert_allclose(alpha.sum() * ratio, gamma @ fractions, rtol=1e-13)
    assert_allclose(
        continuum.spontaneous_recombination, [1.0, 2.0], rtol=0.0, atol=0.0
    )


def test_frozen_beta_initializer_does_not_freeze_nonlinear_solution(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Recompute trapping after the linear seed and recover a scalar root.

    Regime: up=1 s^-1, down=2*beta s^-1, N=1e10 cm^-3 and tau=10*(f0-f1/2).
    An independent scalar balance/escape formula gives the unique root.
    Population rtol=1e-8 follows physical acceptance; beta rtol=1e-8 permits
    propagation through the opacity formula. The frozen answer differs by >0.1.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()
    continuum = ContinuumRateCoefficients(
        np.zeros(2), np.zeros(2), np.ones(2), np.zeros(2), np.zeros(2)
    )
    bound_bound = BoundBoundMatrixRates(
        2,
        np.array([0, 1]),
        np.array([1, 0]),
        np.array([0.0, 2.0]),
        np.zeros(2),
        np.array([1.0, 0.0]),
        np.array([0, 0]),
    )

    def independent_balance(excited_fraction: float) -> float:
        tau = max(10.0 * (1.0 - 1.5 * excited_fraction), 0.0)
        beta = -np.expm1(-tau) / tau if tau > 0.0 else 1.0
        return 1.0 - excited_fraction - 2.0 * beta * excited_fraction

    expected_excited = brentq(independent_balance, 0.0, 2.0 / 3.0, xtol=1e-14)
    expected_tau = 10.0 * (1.0 - 1.5 * expected_excited)
    expected_beta = -np.expm1(-expected_tau) / expected_tau
    initial_beta = -np.expm1(-8.5) / 8.5
    frozen_excited = 1.0 / (1.0 + 2.0 * initial_beta)
    assert abs(frozen_excited - expected_excited) > 0.1
    seed = np.array([0.9, 0.1])
    fractions, _, beta, _ = evaluator._calculate_level_solution(
        0,
        1e9,
        1e4,
        seed,
        continuum,
        evaluator.shell_number_densities[0],
        replace(evaluator.sobolev_inputs[0], tau_coefficient=np.array([1e-9])),
        bound_bound,
    )
    assert_allclose(fractions[1], expected_excited, rtol=1e-8, atol=0.0)
    assert_allclose(beta, [expected_beta], rtol=1e-8, atol=0.0)
    assert abs(beta[0] - initial_beta) > 0.1
    assert (
        abs(independent_balance(fractions[1]))
        <= evaluator.level_residual_tolerance
    )
    assert_allclose(seed, [0.9, 0.1], rtol=0.0, atol=0.0)


@pytest.mark.parametrize(
    "linear_failure", ["singular", "negative", "nonfinite"]
)
def test_unusable_frozen_beta_initializer_keeps_nonlinear_fallback(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    linear_failure: str,
) -> None:
    """An optional linear failure cannot suppress an available physical solve.

    Regime: up/down=2/6 s^-1, exact fractions [.75,.25]. Only the new linear
    boundary is replaced; the real HYBR/log portfolio must solve the atom.
    This is failure-path orchestration, with an independent equilibrium oracle.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()

    def failed_linear_solve(
        *args: object, **kwargs: object
    ) -> npt.NDArray[np.float64]:
        del args, kwargs
        if linear_failure == "singular":
            raise np.linalg.LinAlgError("controlled singular frozen matrix")
        return np.array(
            [1.1, -0.1] if linear_failure == "negative" else [np.nan, np.nan]
        )

    monkeypatch.setattr(np.linalg, "solve", failed_linear_solve)
    continuum = ContinuumRateCoefficients(
        np.zeros(2), np.zeros(2), np.ones(2), np.zeros(2), np.zeros(2)
    )
    bound_bound = BoundBoundMatrixRates(
        2,
        np.array([0, 1]),
        np.array([1, 0]),
        np.zeros(2),
        np.zeros(2),
        np.array([2.0, 6.0]),
        np.full(2, -1),
    )
    with caplog.at_level("DEBUG", logger="tardis.plasma.equilibrium.evaluator"):
        fractions, _, _, _ = evaluator._calculate_level_solution(
            0,
            1e9,
            1e4,
            np.array([0.6, 0.4]),
            continuum,
            evaluator.shell_number_densities[0],
            evaluator.sobolev_inputs[0],
            bound_bound,
        )
    assert_allclose(fractions, [0.75, 0.25], rtol=1e-8, atol=0.0)
    assert "Frozen-beta level initializer" in caplog.text
    assert "retained candidate=HYBR" in caplog.text


def test_disconnected_frozen_beta_matrix_is_skipped(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A disconnected third level does not acquire an invented linear seed.

    Regime: only 0<->1 rates 1/2 s^-1, isolated level 2, no ionization.
    Normalization cannot uniquely close this matrix. The initializer must
    skip it; the existing nonlinear route may retain any physical distribution
    satisfying f0=2*f1. This verifies optional initialization, not uniqueness.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()
    evaluator.level_population_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1), (1, 0, 2)],
        names=evaluator.level_population_index.names,
    )
    seed = np.array([0.1, 0.8, 0.1])
    geometry = replace(
        evaluator.shell_number_densities[0],
        level_number_density=1e10 * seed,
        species_level_positions=np.arange(3),
    )
    before = geometry.level_number_density.copy()
    continuum = ContinuumRateCoefficients(
        np.zeros(3), np.zeros(3), np.ones(3), np.zeros(3), np.zeros(3)
    )
    bound_bound = BoundBoundMatrixRates(
        3,
        np.array([0, 1]),
        np.array([1, 0]),
        np.zeros(2),
        np.zeros(2),
        np.array([1.0, 2.0]),
        np.full(2, -1),
    )
    with caplog.at_level("DEBUG", logger="tardis.plasma.equilibrium.evaluator"):
        fractions, _, _, _ = evaluator._calculate_level_solution(
            0,
            1e9,
            1e4,
            seed,
            continuum,
            geometry,
            evaluator.sobolev_inputs[0],
            bound_bound,
        )
    assert "singular or nonfinite frozen matrix" in caplog.text
    assert_allclose(fractions[0], 2 * fractions[1], rtol=1e-8, atol=0.0)
    assert_allclose(fractions.sum(), 1.0, rtol=0.0, atol=3e-16)
    assert np.all(fractions >= 0.0)
    assert_allclose(geometry.level_number_density, before, rtol=0.0, atol=0.0)
    assert_allclose(seed, [0.1, 0.8, 0.1], rtol=0.0, atol=0.0)


def test_nonlinearly_worse_frozen_seed_is_tried_without_losing_original_best(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A valid frozen seed is tried even when its refreshed beta makes it worse.

    Regime: the nonlinear two-level atom above, original f1=.53. Its frozen
    prediction has f1=.540536 and a larger true relative imbalance. Controlled
    stopped optimizers verify first-trial ownership and original-best retention;
    the physical residual ordering follows independent explicit line fluxes.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()
    original = np.array([0.47, 0.53])
    beta = -np.expm1(-2.05) / 2.05
    predicted_excited = 1.0 / (1.0 + 2.0 * beta)
    root_seeds = []

    def stopped_root(
        _residual: object, seed: npt.NDArray[np.float64], **kwargs: object
    ) -> SimpleNamespace:
        del kwargs
        root_seeds.append(seed.copy())
        return SimpleNamespace(x=original.copy(), success=False)

    def stopped_log(*args: object, **kwargs: object) -> SimpleNamespace:
        del args, kwargs
        return SimpleNamespace(
            x=np.log(original[1:] / original[0]), success=False
        )

    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.root", stopped_root
    )
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.least_squares", stopped_log
    )
    continuum = ContinuumRateCoefficients(
        np.zeros(2), np.zeros(2), np.ones(2), np.zeros(2), np.zeros(2)
    )
    bound_bound = BoundBoundMatrixRates(
        2,
        np.array([0, 1]),
        np.array([1, 0]),
        np.array([0.0, 2.0]),
        np.zeros(2),
        np.array([1.0, 0.0]),
        np.array([0, 0]),
    )
    fractions, _, _, _ = evaluator._calculate_level_solution(
        0,
        1e9,
        1e4,
        original,
        continuum,
        evaluator.shell_number_densities[0],
        replace(evaluator.sobolev_inputs[0], tau_coefficient=np.array([1e-9])),
        bound_bound,
    )
    physical_errors = []
    for excited in (original[1], predicted_excited):
        tau = 10.0 * (1.0 - 1.5 * excited)
        refreshed_beta = -np.expm1(-tau) / tau
        up = 1.0 - excited
        down = 2.0 * refreshed_beta * excited
        physical_errors.append(abs(up - down) / (up + down))
    assert physical_errors[1] > physical_errors[0]
    assert_allclose(root_seeds[0][1], predicted_excited, rtol=1e-13, atol=0.0)
    assert len(root_seeds) == 3  # linear, original, and floored retries survive
    assert_allclose(fractions, original, rtol=1e-14, atol=0.0)
