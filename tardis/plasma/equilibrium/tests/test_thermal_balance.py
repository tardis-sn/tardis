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
    trust-region solve reach the root of a shell with a curved valley and a
    heating row 100 times less sensitive than its electron row, without
    disturbing a well-conditioned shell. This is the inner-shell stall of the
    Type IIP model, where cooling is dominated by collisional pairs close to
    detailed balance.
    Regime: Two shells, solved to the optimizer's own tolerance within its
    evaluation budget.
    Verification: The root is known (T = 1.2, n = T**2). Both residuals of
    both shells must be within 1e-8, the order of the workflow's electron
    limit. Without row scaling the stiff shell stays far from the root at the
    evaluation budget.
    """
    solver = ElectronThermalBalanceSolver(
        CurvedValleyEvaluator(), np.ones(2), pd.DataFrame(), np.ones(2)
    )

    result = solver.solve_multiple_shells(
        np.array([0.8, 0.9, 0.8, 0.9]),
        np.array([0.01, 0.5, 0.01, 0.5]),
        np.array([10.0, 2.0, 10.0, 2.0]),
    )

    np.testing.assert_array_less(
        np.abs(solver.calculate_residual(result.x)),
        1e-8,
        err_msg="a shell did not reach the thermal root",
    )


def test_failed_shell_warning_reports_each_failed_shell(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Name each shell outside its electron and heating limits, with ratios.

    Claim: The closure warning lists exactly the shells with a residual
    outside its limit, with |residual| / limit for both residuals.
    Regime: Three shells against the electron limit 2e-8 and the heating
    limit 2e-7; shell 1 misses the heating limit and shell 2 the electron
    limit.
    Verification: Orchestration check of the warning that reports failed
    shells; the expected ratios follow by dividing the given residuals by the
    limits.
    """
    solver = ElectronThermalBalanceSolver(
        None, np.ones(3), pd.DataFrame(), np.ones(3)
    )
    evaluation = SimpleNamespace(
        electron_residual=pd.Series([1e-8, 1e-8, -6e-8]),
        fractional_heating=pd.Series([2e-8, 4e-7, 0.0]),
    )

    with caplog.at_level(
        "WARNING", logger="tardis.plasma.equilibrium.thermal_balance"
    ):
        solver.warn_failed_shells(evaluation)

    assert "in 2 of 3 shells" in caplog.text
    assert "shell 1: electron 0.5, heating 2" in caplog.text
    assert "shell 2: electron 3, heating 0" in caplog.text
    assert "shell 0:" not in caplog.text
