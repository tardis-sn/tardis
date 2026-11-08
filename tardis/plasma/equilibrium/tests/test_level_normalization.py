"""Verify population normalization independently of dimensional rate sources."""

from dataclasses import replace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from scipy.optimize import root

from tardis.plasma.equilibrium.evaluator import (
    calculate_nlte_level_population_residual,
)
from tardis.plasma.equilibrium.inputs import (
    BoundBoundMatrixRates,
    LevelEquationRates,
    ShellNumberDensity,
    SobolevInputs,
)


@pytest.fixture
def two_level_normalization_inputs() -> tuple[
    ShellNumberDensity, SobolevInputs, BoundBoundMatrixRates
]:
    """Supply two bound levels with a fixed 6 s^-1 downward transition."""
    empty_integer = np.array([], dtype=np.int64)
    empty_float = np.array([], dtype=np.float64)
    empty_boolean = np.array([], dtype=np.bool_)
    shell_number_density = ShellNumberDensity(
        1.0, np.array([0.5, 0.5]), np.array([0, 1])
    )
    sobolev = SobolevInputs(
        empty_integer,
        empty_integer,
        empty_float,
        empty_float,
        empty_boolean,
        empty_boolean,
        empty_float,
        empty_integer,
        pd.Index([]),
    )
    bound_bound_rates = BoundBoundMatrixRates(
        2,
        np.array([1]),
        np.array([0]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([6.0]),
        np.array([-1]),
    )
    return shell_number_density, sobolev, bound_bound_rates


@pytest.mark.parametrize("case_b", [False, True])
@pytest.mark.parametrize("rate_scale", [1.0e-6, 1.0, 1.0e6])
def test_normalization_row_excludes_recombination_source(
    two_level_normalization_inputs: tuple[
        ShellNumberDensity, SobolevInputs, BoundBoundMatrixRates
    ],
    case_b: bool,
    rate_scale: float,
) -> None:
    """Claim: the constraint is sum(b)-1, not a dimensionful kinetic equation.

    Regime: two levels, ionization [2, 4] and recombination [3, 5] s^-1,
    scaled together over twelve orders of magnitude. An off-normalization
    candidate must retain its 0.1 defect under both recombination policies.
    Verification is the exact sum constraint; atol=2e-16 allows roundoff in
    summing two fractions and subtracting one, independently of rate strength.
    """
    shell_number_density, sobolev, bound_bound_rates = (
        two_level_normalization_inputs
    )
    ionization = np.array([2.0, 4.0]) * rate_scale
    level_rates = LevelEquationRates(
        ionization,
        np.array([3.0, 5.0]) * rate_scale,
        -np.diag(ionization),
    )
    residual, _, ratio = calculate_nlte_level_population_residual(
        np.array([0.6, 0.5]),
        level_rates,
        None,
        None,
        None,
        (2, 0),
        shell_number_density,
        sobolev,
        bound_bound_rates=bound_bound_rates,
        case_b=case_b,
    )

    np.testing.assert_allclose(residual[0], 0.1, rtol=0.0, atol=2.0e-16)
    # Case A must still include direct ground recombination in ion balance.
    np.testing.assert_allclose(
        ratio, 3.2 / (5.0 if case_b else 8.0), rtol=2.0e-15
    )


@pytest.mark.parametrize("case_b", [False, True])
@pytest.mark.parametrize("rate_scale", [1.0e-6, 1.0, 1.0e6])
def test_reduced_two_level_solution_preserves_normalization(
    two_level_normalization_inputs: tuple[
        ShellNumberDensity, SobolevInputs, BoundBoundMatrixRates
    ],
    case_b: bool,
    rate_scale: float,
) -> None:
    """Claim: normalized equilibrium fractions do not depend on the rate unit.

    Regime: downward rate 6, ionization [2, 4], recombination [3, 5] s^-1;
    Case B excludes the ground recombination 3. Analytic excited-state balance
    gives b0=6*b1 in Case A and b0=3*b1 in Case B. Thus b=[6/7,1/7] or
    [3/4,1/4]. The independent oracle combines those balances with sum(b)=1.
    rtol=1e-12 covers floating-point root solution of this linear system.
    """
    shell_number_density, sobolev, bound_bound_rates = (
        two_level_normalization_inputs
    )
    ionization = np.array([2.0, 4.0]) * rate_scale
    level_rates = LevelEquationRates(
        ionization,
        np.array([3.0, 5.0]) * rate_scale,
        -np.diag(ionization),
    )
    bound_bound_rates = replace(
        bound_bound_rates, fixed_rate=bound_bound_rates.fixed_rate * rate_scale
    )

    def calculate_residual(
        fractions: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return calculate_nlte_level_population_residual(
            fractions,
            level_rates,
            None,
            None,
            None,
            (2, 0),
            shell_number_density,
            sobolev,
            bound_bound_rates=bound_bound_rates,
            case_b=case_b,
        )[0]

    solution = root(calculate_residual, np.array([0.5, 0.5]))
    expected_fractions = [0.75, 0.25] if case_b else [6.0 / 7.0, 1.0 / 7.0]
    np.testing.assert_allclose(
        solution.x, expected_fractions, rtol=1.0e-12, atol=0.0
    )
    np.testing.assert_allclose(solution.x.sum(), 1.0, rtol=0.0, atol=2.0e-15)
