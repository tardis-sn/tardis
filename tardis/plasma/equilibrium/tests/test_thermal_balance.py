from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest

from tardis.plasma.equilibrium.thermal_balance import (
    ElectronThermalBalanceSolver,
    calculate_batch_standard_error,
)
from tardis.plasma.exceptions import PlasmaIonizationError


def test_batch_standard_error_matches_pooled_monte_carlo_noise() -> None:
    """Recover the noise of a pooled estimate from unequal packet batches.

    Claim: Batches holding n_b packets of a per-packet scatter sigma give the
    standard error sigma / sqrt(N) of the packet-weighted mean of all N
    packets.
    Regime: 16 batches of unequal size, as uneven Monte Carlo threads give.
    Verification: Batch means are drawn with variance sigma**2 / n_b, so the
    exact answer is known. Averaging the variance estimate over many trials
    removes its 16-batch sampling scatter.
    """
    rng = np.random.default_rng(23111963)
    per_packet_sigma = 3.0
    packet_counts = rng.integers(500, 1500, size=16)
    expected_variance = per_packet_sigma**2 / packet_counts.sum()

    trials = 4000
    batch_means = 10.0 + rng.normal(
        size=(trials, len(packet_counts))
    ) * per_packet_sigma / np.sqrt(packet_counts)
    estimated_variance = np.mean(
        [
            calculate_batch_standard_error(trial_means, packet_counts) ** 2
            for trial_means in batch_means
        ]
    )

    # One 16-batch variance estimate scatters by sqrt(2 / 15) ~ 37 %; the
    # mean of 4000 scatters by ~0.6 %, so 3 % is about five sigma.
    np.testing.assert_allclose(estimated_variance, expected_variance, rtol=0.03)


class CurvedValleyEvaluator:
    """Two shells with a curved, stiff electron-temperature valley.

    Shell 0 mimics an inner shell: its heating residual responds 100 times
    less to the temperature than its electron residual does to the density.
    Shell 1 is well conditioned. The density and temperature are the solver's
    candidate fractions, since the maximum density and radiation temperature
    are 1.
    """

    heating_sensitivity = np.array([1e-2, 1.0])
    root_temperature = 1.2

    def evaluate(
        self,
        electron_density: npt.NDArray[np.float64],
        electron_temperature: npt.NDArray[np.float64],
        level_seed: pd.DataFrame,
    ) -> SimpleNamespace:
        return SimpleNamespace(
            electron_residual=pd.Series(
                10.0 * (electron_density - electron_temperature**2)
            ),
            fractional_heating=pd.Series(
                self.heating_sensitivity
                * (self.root_temperature - electron_temperature)
            ),
        )


def test_row_scaling_closes_stiff_valley_shell() -> None:
    """Close a shell whose heating row is much weaker than its electron row.

    Claim: Dividing each residual row by its Jacobian diagonal lets the joint
    trust-region solve close a shell with a curved valley and a heating row 100
    times less sensitive than its electron row, without disturbing a
    well-conditioned shell. This is the inner-shell stall of the Type IIP
    model, where cooling is dominated by collisional pairs in detailed balance.
    Regime: Two shells; limits 1e-4 on the electron residual and 1e-3 of each
    heating sensitivity on the heating residual.
    Verification: The root is known (T = 1.2, n = T**2). Without row scaling
    the stiff shell stays about 50 times over its heating limit at the
    evaluation budget.
    """
    # The noise measurement is not exercised, so the batches are unused.
    solver = ElectronThermalBalanceSolver(
        CurvedValleyEvaluator(),
        np.ones(2),
        pd.DataFrame(),
        np.ones(2),
        batch_evaluators=lambda: (),
        batch_packet_counts=np.ones(2, dtype=np.int64),
    )
    heating_limits = 1e-3 * CurvedValleyEvaluator.heating_sensitivity
    residual_limits = np.array(
        [1e-4, heating_limits[0], 1e-4, heating_limits[1]]
    )

    candidate = solver.solve_multiple_shells(
        np.array([0.8, 0.9, 0.8, 0.9]),
        np.array([0.01, 0.5, 0.01, 0.5]),
        np.array([10.0, 2.0, 10.0, 2.0]),
        residual_limits,
    )

    closed_shells = solver.find_converged_shells(
        solver.calculate_residual(candidate), residual_limits
    )
    assert closed_shells.all(), "a shell did not reach its limits"


def test_open_shell_warning_reports_each_open_shell(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Warn about each shell outside its Monte Carlo limits with its ratio.

    Claim: The closure warning lists exactly the shells with a residual
    outside its limit, with |residual| / limit for both residuals.
    Regime: Three shells with limits of 0.1; shell 1 has an open heating
    residual and shell 2 an open electron residual.
    Verification: Orchestration check of the warning that reports open
    shells; the expected ratios follow by dividing the given residuals by the
    limit.
    """
    solver = ElectronThermalBalanceSolver(
        None,
        np.ones(3),
        pd.DataFrame(),
        np.ones(3),
        batch_evaluators=lambda: (),
        batch_packet_counts=np.ones(2, dtype=np.int64),
    )
    solver.outer_residual_limits = np.full(6, 0.1)
    evaluation = SimpleNamespace(
        electron_residual=pd.Series([0.05, 0.05, -0.3]),
        fractional_heating=pd.Series([0.01, 0.2, 0.0]),
    )

    with caplog.at_level(
        "WARNING", logger="tardis.plasma.equilibrium.thermal_balance"
    ):
        solver.warn_failed_shells(evaluation)

    assert "in 2 of 3 shells" in caplog.text
    assert "shell 1: electron 0.5, heating 2" in caplog.text
    assert "shell 2: electron 3, heating 0" in caplog.text
    assert "shell 0:" not in caplog.text


class BatchEvaluator:
    """Return fixed residuals for one packet batch, or fail its ion solve."""

    def __init__(self, electron_residual: float | None) -> None:
        self.electron_residual = electron_residual

    def evaluate(self, *_args: object) -> SimpleNamespace:
        if self.electron_residual is None:
            raise PlasmaIonizationError("batch ionization failed")
        return SimpleNamespace(
            electron_residual=pd.Series([self.electron_residual]),
            fractional_heating=pd.Series([0.0]),
        )


def test_failed_batch_is_dropped_from_noise_estimate() -> None:
    """Measure the noise from the batches whose ionization balance solves.

    Claim: A packet batch whose ionization balance fails does not abort the
    noise estimate; the standard error comes from the remaining batches and
    their packet counts.
    Regime: One shell and four equal batches, the second of which fails.
    Verification: The three remaining batch residuals 0.1, 0.2 and 0.6 have a
    pooled mean of 0.3 and a standard error of
    sqrt(((0.2**2 + 0.1**2 + 0.3**2) / 2) / 3) = sqrt(0.07 / 3).
    """
    solver = ElectronThermalBalanceSolver(
        None,
        np.ones(1),
        pd.DataFrame(),
        np.ones(1),
        batch_evaluators=lambda: (
            BatchEvaluator(value) for value in (0.1, None, 0.2, 0.6)
        ),
        batch_packet_counts=np.full(4, 100, dtype=np.int64),
    )

    noise = solver.calculate_residual_noise(np.ones(2))

    np.testing.assert_allclose(noise, [np.sqrt(0.07 / 3), 0.0], rtol=1e-12)
