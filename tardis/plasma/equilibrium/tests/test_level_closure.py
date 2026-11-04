"""Physical closure requirements independent of nonlinear solver status."""

from collections.abc import Callable
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from astropy import units as u
from numpy.testing import assert_allclose
from scipy.optimize import brentq

from tardis import constants as const
from tardis.plasma.equilibrium.evaluator import (
    PlasmaEquilibriumEvaluator,
    calculate_nlte_level_population_residual,
)
from tardis.plasma.equilibrium.inputs import (
    BoundBoundMatrixRates,
    ContinuumRateCoefficients,
    LevelEquationRates,
    ShellNumberDensity,
)
from tardis.plasma.equilibrium.rate_matrix import RateMatrix
from tardis.plasma.equilibrium.rates import RadiativeRatesSolver
from tardis.plasma.equilibrium.solver import PlasmaEquilibriumSolver
from tardis.plasma.equilibrium.tests.test_evaluator import (
    toy_evaluator as toy_evaluator,  # noqa: PLC0414 - re-export shared pytest fixture
)


@pytest.fixture
def coupled_helium_evaluator(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> PlasmaEquilibriumEvaluator:
    """Construct a synthetic two-level He I atom with He II/III charge balance.

    Photoionization is [1,3] s^-1; recombination is [1,2]e-10 cm^3/s.
    He III/He II = 1e9/n_e, and N_He=1e10 cm^-3. Line rates are beta and
    2beta with tau=1e-9*(n0-n1/2). Coefficient and ion-solver boundaries are
    analytic so this exercises real reduced level, opacity, and outer solves
    without requiring an atomic database.
    """
    level_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0), (2, 0, 1)],
        names=toy_evaluator.level_population_index.names,
    )
    all_levels = pd.MultiIndex.from_tuples(
        [(2, 0, 0), (2, 0, 1), (2, 1, 0), (2, 2, 0)],
        names=level_index.names,
    )
    levels = pd.DataFrame(
        {"energy": [0.0, 1e-12, 0.0, 0.0], "g": [2.0, 4.0, 1.0, 1.0]},
        index=all_levels,
    )
    ion_index = pd.MultiIndex.from_tuples(
        [(2, 0), (2, 1), (2, 2)], names=["atomic_number", "ion_number"]
    )
    line_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0, 1)], names=toy_evaluator.j_blues.index.names
    )
    cross_sections = toy_evaluator.photoionization_cross_sections.copy()
    cross_sections.index = level_index
    evaluator = PlasmaEquilibriumEvaluator(
        cross_sections,
        pd.Series([0, 1], index=level_index),
        toy_evaluator.estimators_continuum,
        2e5 * u.s,
        3e30 * u.cm**3,
        levels,
        pd.Series([2e-11, 4e-11], index=ion_index[1:]),
        RateMatrix(
            RadiativeRatesSolver(
                pd.DataFrame(
                    {"A_ul": [2.0], "B_ul": [0.0], "B_lu": [1.0], "nu": [1.0]},
                    index=line_index,
                )
            ),
            toy_evaluator.rate_matrix_solver.electron_rate_solver,
            levels,
        ),
        pd.DataFrame([1.0], index=line_index),
        (ShellNumberDensity(1e10, np.array([1e10, 0.0]), np.array([0, 1])),),
        (
            replace(
                toy_evaluator.sobolev_inputs[0],
                line_index=line_index,
                tau_coefficient=np.array([1e-9]),
            ),
        ),
        level_index,
        (2, 0),
        pd.DataFrame([[1e10]], index=[2], columns=[0]),
        np.array([2e10]),
        case_b_species=(),
    )

    class AnalyticCoefficients:
        def solve_coefficients(
            self,
            temperature: u.Quantity,
        ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            kelvin = temperature.to_value(u.K)[0]
            saha_inverse = (
                np.array([1.0, 2.0])
                * (
                    const.h.cgs.value**2
                    / (
                        2
                        * np.pi
                        * const.m_e.cgs.value
                        * const.k_B.cgs.value
                        * kelvin
                    )
                )
                ** 1.5
                * np.exp(
                    np.array([2e-11, 1.9e-11]) / (const.k_B.cgs.value * kelvin)
                )
            )
            return (
                pd.DataFrame([1.0, 3.0], index=level_index),
                pd.DataFrame([0.0, 0.0], index=level_index),
                pd.DataFrame(
                    np.array([1e-10, 2e-10]) / saha_inverse, index=level_index
                ),
            )

    class ZeroCollisionalIonization:
        def solve_coefficients(self, temperature: u.Quantity) -> pd.DataFrame:
            del temperature
            return pd.DataFrame([0.0, 0.0], index=level_index)

    class AnalyticIonSolver:
        def solve(
            self,
            estimated_level_population: pd.DataFrame,
            **kwargs: object,
        ) -> tuple[pd.DataFrame, pd.Series]:
            del kwargs
            estimated = estimated_level_population.iloc[:, 0].to_numpy()
            fractions = estimated / estimated.sum()
            ionization = np.dot([1.0, 3.0], fractions)

            def charge_error(density: float) -> float:
                ratio = ionization / (3e-10 * density)
                third_ratio = 1e9 / density
                neutral = 1e10 / (1 + ratio * (1 + third_ratio))
                return neutral * ratio * (1 + 2 * third_ratio) - density

            density = brentq(charge_error, 1.0, 2e10, xtol=1e-5, rtol=1e-14)
            ratio = ionization / (3e-10 * density)
            neutral = 1e10 / (1 + ratio * (1 + 1e9 / density))
            return (
                pd.DataFrame(
                    [neutral, neutral * ratio, neutral * ratio * 1e9 / density],
                    index=ion_index,
                    columns=[0],
                ),
                pd.Series([density], index=[0]),
            )

    evaluator.estimated_photoionization_rate_solver = AnalyticCoefficients()
    evaluator.collisional_ionization_rate_solver = ZeroCollisionalIonization()
    evaluator.ion_population_solver = AnalyticIonSolver()
    return evaluator


def test_refinement_excludes_third_stage_from_opacity_budget_and_copies_snapshot(
    coupled_helium_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """He III nuclei cannot enter the reduced He I/II opacity budget.

    Nuclear conservation gives N_pair=N_He-n_HeIII independently of the
    reduced solver. Refinement must copy the supplied authoritative snapshot,
    not mutate it or retain all helium as the adjacent-pair budget.
    """
    evaluator = coupled_helium_evaluator
    seed = pd.DataFrame(
        [[0.8], [0.2]], index=evaluator.level_population_index, columns=[0]
    )
    previous = evaluator.evaluate(np.array([7e9]), np.array([9000.0]), seed)
    levels_before = previous.absolute_level_population.copy(deep=True)
    ions_before = previous.ion_population.copy(deep=True)
    seed_before = seed.copy(deep=True)
    evaluator.refine_terminal(np.array([7e9]), np.array([9000.0]), previous)
    assert ions_before.loc[(2, 2), 0] > 1e8
    assert_allclose(
        evaluator.shell_number_densities[0].species_number_density,
        1e10 - ions_before.loc[(2, 2), 0],
        rtol=1e-14,
    )
    assert_allclose(
        evaluator.shell_number_densities[0].level_number_density,
        levels_before.iloc[:, 0],
        rtol=0.0,
        atol=0.0,
    )
    pd.testing.assert_frame_equal(
        previous.absolute_level_population, levels_before
    )
    pd.testing.assert_frame_equal(previous.ion_population, ions_before)
    pd.testing.assert_frame_equal(seed, seed_before)


def test_terminal_refinement_closes_real_three_stage_opacity_coupling(
    coupled_helium_evaluator: PlasmaEquilibriumEvaluator,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Refine real beta/ion coupling and recover an independent scalar root.

    The initial eliminated continuum wrongly assigns He III's nuclei to the
    He I/II pair. The final coupled solution is independently obtained by
    scalar charge balance and the analytic excited equation
    beta*(b0-2*b1)+(2/3)*b0-b1=0. Population tolerance 1e-7 permits the
    production 1e-8 kinetic target and nonlinear conditioning; all returned
    physical residuals must still satisfy their unchanged solver tolerances.
    """
    evaluator = coupled_helium_evaluator
    seed = pd.DataFrame(
        [[0.8], [0.2]], index=evaluator.level_population_index, columns=[0]
    )
    initial = evaluator.evaluate(np.array([7e9]), np.array([9000.0]), seed)
    assert np.max(abs(initial.relative_level_residual.to_numpy())) > 1e-4

    def excited_balance(excited: float) -> float:
        ionization = 1 + 2 * excited

        def charge_error(density: float) -> float:
            ratio = ionization / (3e-10 * density)
            neutral = 1e10 / (1 + ratio * (1 + 1e9 / density))
            return neutral * ratio * (1 + 2e9 / density) - density

        density = brentq(charge_error, 1.0, 2e10, rtol=1e-14)
        ratio = ionization / (3e-10 * density)
        neutral = 1e10 / (1 + ratio * (1 + 1e9 / density))
        tau = 1e-9 * neutral * (1 - 1.5 * excited)
        beta = -np.expm1(-tau) / tau
        return beta * (1 - 3 * excited) + 2 / 3 - 5 / 3 * excited

    expected_excited = brentq(excited_balance, 0.0, 0.6, xtol=1e-14)
    caplog.set_level("INFO", logger="tardis.plasma.equilibrium.solver")
    solution = PlasmaEquilibriumSolver(fixed_temperature_ratio=0.9).solve(
        evaluator,
        [7e9],
        [9000.0],
        seed,
        [2e10],
        [1e4],
    )
    final = solution.evaluation
    assert_allclose(
        final.normalized_population.iloc[1, 0], expected_excited, rtol=1e-7
    )
    assert (
        np.max(abs(final.relative_level_residual.to_numpy()))
        <= evaluator.level_residual_tolerance
    )
    assert_allclose(final.ion_population.sum(), 1e10, rtol=1e-14)
    charge = final.ion_population.multiply([0, 1, 2], axis=0).sum()
    assert_allclose(charge, solution.electron_density, rtol=1e-7)
    assert "Starting terminal level refinement" in caplog.text
    # The context now contains the adjacent pair, not all helium nuclei.
    assert evaluator.shell_number_densities[0].species_number_density < 1e10


@pytest.mark.parametrize("rate_scale", [1e-24, 1.0, 1e24])
def test_reduced_levels_are_invariant_to_uniform_rate_scale(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    rate_scale: float,
) -> None:
    """Resolve the same analytic equilibrium at arbitrarily small rates.

    Two optically thin levels have up/down rates 1/2, ionization rates
    [1,3], and recombination rates [1,2] per upper ion. Eliminating the
    continuum gives (5/3)b0-3b1=0, hence b=[9/14,5/14]. The configured
    1e-8 relative kinetic target determines population tolerance; an
    independent net/gross calculation tests physical rather than optimizer
    closure. Scaling every rate cannot change either criterion.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()
    sobolev = replace(evaluator.sobolev_inputs[0], tau_coefficient=np.zeros(1))
    bound_bound = BoundBoundMatrixRates(
        2,
        np.array([0, 1]),
        np.array([1, 0]),
        rate_scale * np.array([1.0, 2.0]),
        np.zeros(2),
        np.zeros(2),
        np.array([0, 0]),
    )
    continuum = ContinuumRateCoefficients(
        rate_scale * np.array([1.0, 3.0]),
        np.zeros(2),
        rate_scale * np.array([1.0, 2.0]) / 1e9,
        np.zeros(2),
        np.zeros(2),
    )
    fractions, ratio, _, _ = evaluator._calculate_level_solution(
        0,
        1e9,
        1e4,
        np.array([0.9, 0.1]),
        continuum,
        evaluator.shell_number_densities[0],
        sobolev,
        bound_bound,
    )
    assert_allclose(fractions, [9 / 14, 5 / 14], rtol=1e-8, atol=0.0)
    matrix = np.array([[-2.0, 2.0], [1.0, -5.0]])
    source = np.array([1.0, 2.0]) * ratio
    net = matrix @ fractions + source
    gross = abs(matrix) @ fractions + source
    assert np.max(abs(net) / gross) <= evaluator.level_residual_tolerance


def test_best_effort_candidate_uses_relative_physical_rank(
    toy_evaluator: PlasmaEquilibriumEvaluator,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Retain the best relative candidate even if HYBR finally returns worse.

    Optically thin up/down rates 1/100 give net excited flow b0-100b1.
    Candidate [0.98,0.02] has larger absolute error (1.02) but smaller
    relative error (1.02/2.98) than [0.999,0.001] (0.899/1.099).
    Controlled SciPy stops exercise best-effort and callback wiring only;
    the independent analytic ordering is the numerical oracle.
    """
    better = np.array([0.98, 0.02])
    worse = np.array([0.999, 0.001])

    def stopped_root(
        objective: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
        seed: npt.NDArray[np.float64],
        **kwargs: object,
    ) -> SimpleNamespace:
        del seed, kwargs
        objective(better)
        objective(worse)
        return SimpleNamespace(x=worse, success=True)

    def stopped_log(
        objective: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
        seed: npt.NDArray[np.float64],
        **kwargs: object,
    ) -> SimpleNamespace:
        del objective, seed, kwargs
        return SimpleNamespace(x=np.log(worse[1:] / worse[0]), success=False)

    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.root", stopped_root
    )
    monkeypatch.setattr(
        "tardis.plasma.equilibrium.evaluator.least_squares", stopped_log
    )
    evaluator = toy_evaluator
    bound_bound = BoundBoundMatrixRates(
        2,
        np.array([0, 1]),
        np.array([1, 0]),
        np.array([1.0, 100.0]),
        np.zeros(2),
        np.zeros(2),
        np.array([0, 0]),
    )
    continuum = ContinuumRateCoefficients(
        np.zeros(2),
        np.zeros(2),
        np.ones(2),
        np.zeros(2),
        np.zeros(2),
    )
    fractions, _, _, _ = evaluator._calculate_level_solution(
        0,
        1e9,
        1e4,
        np.array([0.5, 0.5]),
        continuum,
        evaluator.shell_number_densities[0],
        replace(evaluator.sobolev_inputs[0], tau_coefficient=np.zeros(1)),
        bound_bound,
    )
    assert_allclose(fractions, better, rtol=0.0, atol=0.0)
    # The best physical trial is deliberately still open; intermediate
    # best-effort must not become an unconditional strict inner rejection.
    assert abs(fractions[0] - 100 * fractions[1]) > 1.0


def test_equilibrium_terminal_validation_rejects_unclosed_levels() -> None:
    """Exact charge balance cannot compensate for a 20% kinetic imbalance.

    This is acceptance wiring, not a kinetic solver verification. The reported
    relative residual is independently prescribed while charge is exactly zero.
    """

    class UnclosedEvaluator:
        def evaluate(
            self,
            electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            return SimpleNamespace(
                electron_residual=pd.Series([0.0]),
                normalized_population=level_seed,
                relative_level_residual=pd.DataFrame([[0.0], [0.2]]),
            )

    with pytest.raises(RuntimeError, match=r"level.*residual"):
        PlasmaEquilibriumSolver(fixed_temperature_ratio=0.9).solve(
            UnclosedEvaluator(),
            [1e9],
            [9000.0],
            pd.DataFrame([[0.5], [0.5]]),
            [1e10],
            [1e4],
        )


def test_final_level_closure_uses_absolute_populations_and_continuum_ratio(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Rebuild the physical final rates, not trial fractions or trial beta.

    Two levels with b=[0.8,0.2], n_ion/n_neutral=0.25, ionization 1 s^-1,
    Case B recombination 1 s^-1, and radiative up/down rates beta, 2beta.
    The final optical depth is exactly 7, and the excited balance is
    0.4*beta+0.05 with gross flow 1.2*beta+0.45. The trial seed [0.5,0.5]
    and its eliminated continuum ratio differ deliberately. Atol=0 ensures
    the positive imbalance cannot be hidden; rtol=1e-12 allows roundoff.
    """
    evaluator = toy_evaluator
    evaluator.sobolev_inputs = (
        replace(evaluator.sobolev_inputs[0], tau_coefficient=np.array([1e-9])),
    )
    levels = pd.DataFrame(
        [[8e9], [2e9]], index=evaluator.level_population_index, columns=[0]
    )
    ions = pd.DataFrame(
        [[1e10], [2.5e9]],
        index=pd.MultiIndex.from_tuples([(1, 0), (1, 1)]),
        columns=[0],
    )
    continuum = ContinuumRateCoefficients(
        np.ones(2), np.zeros(2), np.full(2, 1e-9), np.zeros(2), np.zeros(2)
    )
    state = SimpleNamespace(
        fractions=np.array([[0.5, 0.5]]),
        species_index=evaluator.level_population_index,
        continuum_rate_coefficients=(continuum,),
        bound_bound_rates=(None,),
    )
    residual = evaluator._calculate_final_residual(
        pd.Series([1e9]),
        np.array([1e4]),
        levels,
        state,
        state,
        normalize_residual=True,
        ion_population=ions,
    )
    beta = -np.expm1(-7.0) / 7.0
    denominator = (
        1.2 * beta
        + 0.45
        + 64
        * np.finfo(float).eps
        / evaluator.level_residual_tolerance
        * (1.25 + 2.4 * beta)
    )
    ground_denominator = (
        0.8
        + 1.2 * beta
        + 64
        * np.finfo(float).eps
        / evaluator.level_residual_tolerance
        * (1.25 + 2.4 * beta)
    )
    assert_allclose(
        residual.to_numpy()[:, 0],
        [
            (-0.8 - 0.4 * beta) / ground_denominator,
            (0.4 * beta + 0.05) / denominator,
        ],
        rtol=1e-12,
        atol=0.0,
    )


def test_rate_normalization_has_zero_turnover_limit(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Absent transitions imply zero kinetic imbalance, not 0/0.

    Verification is the exact no-flow limiting case: zero ionization and
    bound-bound rates produce zero continuum ratio and source. Recombination
    coefficients are positive but no upper ions are populated. Any normalized
    bound distribution conserves population. Absolute tolerance is zero.
    """
    empty_integer = np.array([], dtype=np.int64)
    empty_float = np.array([], dtype=np.float64)
    rates = BoundBoundMatrixRates(
        2,
        empty_integer,
        empty_integer,
        empty_float,
        empty_float,
        empty_float,
        empty_integer,
    )
    residual, _, ratio = calculate_nlte_level_population_residual(
        np.array([0.6, 0.4]),
        LevelEquationRates(np.zeros(2), np.ones(2), np.zeros((2, 2))),
        toy_evaluator.rate_matrix_solver,
        None,
        None,
        (1, 0),
        toy_evaluator.shell_number_densities[0],
        toy_evaluator.sobolev_inputs[0],
        bound_bound_rates=rates,
        normalize_residual=True,
    )
    assert_allclose(residual, 0.0, atol=0.0)
    assert ratio == 0.0


def test_final_closure_retains_ground_kinetics_when_excited_levels_close(
    toy_evaluator: PlasmaEquilibriumEvaluator,
) -> None:
    """Ground-only continuum imbalance must not disappear into normalization.

    With no upward pumping, all neutral atoms are in the ground state and
    the excited equation is exactly closed. Ground ionization is 1 s^-1
    while recombination supplies only 0.25 s^-1 per neutral atom. Its net
    -0.75 s^-1 and gross 1.25 s^-1 must remain in final kinetic validation.
    Verification uses these analytic flows with roundoff-only rtol=1e-14.
    """
    evaluator = toy_evaluator
    evaluator.case_b_species = frozenset()
    evaluator.j_blues.iloc[:, :] = 0.0
    levels = pd.DataFrame(
        [[1e10], [0.0]], index=evaluator.level_population_index, columns=[0]
    )
    ions = pd.DataFrame(
        [[1e10], [2.5e9]],
        index=pd.MultiIndex.from_tuples([(1, 0), (1, 1)]),
        columns=[0],
    )
    continuum = ContinuumRateCoefficients(
        np.array([1.0, 0.0]),
        np.zeros(2),
        np.array([1e-9, 0.0]),
        np.zeros(2),
        np.zeros(2),
    )
    state = SimpleNamespace(
        fractions=np.array([[1.0, 0.0]]),
        species_index=evaluator.level_population_index,
        continuum_rate_coefficients=(continuum,),
        bound_bound_rates=(None,),
    )
    actual = evaluator._calculate_final_residual(
        pd.Series([1e9]),
        np.array([1e4]),
        levels,
        state,
        state,
        normalize_residual=True,
        ion_population=ions,
    )
    scale = 1.25 * (
        1 + 64 * np.finfo(float).eps / evaluator.level_residual_tolerance
    )
    assert_allclose(
        actual.to_numpy()[:, 0], [-0.75 / scale, 0.0], rtol=1e-14, atol=0.0
    )
