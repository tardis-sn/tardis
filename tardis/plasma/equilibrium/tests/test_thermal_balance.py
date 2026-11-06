from types import SimpleNamespace

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest

from tardis.plasma.equilibrium.thermal_balance import (
    ElectronThermalBalanceSolver,
)


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
    solver = ElectronThermalBalanceSolver(
        CurvedValleyEvaluator(), np.ones(2), pd.DataFrame(), np.ones(2)
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
    """Name each shell outside its electron and heating limits, with ratios.

    Claim: The closure warning lists exactly the shells with a residual
    outside its limit, with |residual| / limit for both residuals.
    Regime: Three shells with limits of 0.1; shell 1 has an open heating
    residual and shell 2 an open electron residual.
    Verification: Orchestration check of the warning that reports open
    shells; the expected ratios follow by dividing the given residuals by the
    limit.
    """
    solver = ElectronThermalBalanceSolver(
        None, np.ones(3), pd.DataFrame(), np.ones(3)
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
