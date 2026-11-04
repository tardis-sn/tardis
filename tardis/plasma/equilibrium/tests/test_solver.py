import re
from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest

from tardis.plasma.equilibrium.solver import PlasmaEquilibriumSolver


class AnalyticEvaluator:
    """Evaluator with a prescribed density and temperature root."""

    def __init__(self) -> None:
        self.calls = 0

    def evaluate(
        self,
        electron_density: np.ndarray,
        electron_temperature: np.ndarray,
        level_seed: pd.DataFrame,
    ) -> SimpleNamespace:
        self.calls += 1
        columns = level_seed.columns
        normalized_population = level_seed / level_seed.sum(axis=0)
        return SimpleNamespace(
            electron_residual=pd.Series(
                electron_density / 2.0e9 - 1.0, index=columns
            ),
            fractional_heating=pd.Series(
                electron_temperature / 8000.0 - 1.0, index=columns
            ),
            normalized_population=normalized_population,
        )


def test_equilibrium_solver_finds_normalized_density_temperature_root() -> None:
    """Solve electron and thermal residuals and rebuild with accepted levels."""
    evaluator = AnalyticEvaluator()
    level_seed = pd.DataFrame([[0.8], [0.2]], columns=[0])

    solution = PlasmaEquilibriumSolver().solve(
        evaluator,
        np.array([1.0e9]),
        np.array([6000.0]),
        level_seed,
        np.array([1.0e10]),
        np.array([1.0e4]),
    )

    npt.assert_allclose(solution.electron_density, [2.0e9], rtol=1e-10)
    npt.assert_allclose(solution.electron_temperature, [8000.0], rtol=1e-10)
    npt.assert_allclose(solution.evaluation.electron_residual, 0.0, atol=1e-10)
    npt.assert_allclose(solution.evaluation.fractional_heating, 0.0, atol=1e-10)
    assert evaluator.calls > 2


@pytest.mark.parametrize("thermal_residual", [None, np.nan, 1.0e6])
@pytest.mark.parametrize("shell_count", [1, 2])
def test_equilibrium_solver_prescribes_temperature_and_solves_only_charge(
    thermal_residual: float | None,
    shell_count: int,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Prescribe 0.9 T_rad while retaining the exact electron-density root.

    Claim: unused or nonfinite thermal residuals cannot affect a prescribed
    temperature. Regime: two independent linear charge roots and a thermal
    root inconsistent with 0.9 T_rad. Verification: every evaluation, including
    terminal validation, uses the prescribed temperatures; the charge root is
    2e9 cm^-3 to the existing analytic-root tolerance (1e-10 relative).
    """
    temperatures = []

    class FixedTemperatureEvaluator(AnalyticEvaluator):
        def evaluate(
            self,
            electron_density: np.ndarray,
            electron_temperature: np.ndarray,
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            temperatures.append(electron_temperature.copy())
            evaluation = super().evaluate(
                electron_density, electron_temperature, level_seed
            )
            evaluation.fractional_heating = (
                None
                if thermal_residual is None
                else pd.Series(thermal_residual, index=level_seed.columns)
            )
            return evaluation

    caplog.set_level("INFO", logger="tardis.plasma.equilibrium.solver")
    solver = PlasmaEquilibriumSolver(
        fixed_temperature_ratio=0.9, minimum_electron_temperature=9500.0
    )
    shell_columns = [17, 54][:shell_count]
    solution = solver.solve(
        FixedTemperatureEvaluator(),
        np.array([1.0e9, 3.0e9])[:shell_count],
        np.array([6000.0, 7000.0])[:shell_count],
        pd.DataFrame(
            np.tile([[0.8], [0.2]], shell_count), columns=shell_columns
        ),
        np.full(shell_count, 1.0e10),
        np.array([1.0e4, 1.2e4])[:shell_count],
    )

    npt.assert_allclose(solution.electron_density, 2.0e9, rtol=1e-10)
    expected_temperature = np.array([9000.0, 10800.0])[:shell_count]
    npt.assert_array_equal(solution.electron_temperature, expected_temperature)
    for temperature in temperatures:
        npt.assert_array_equal(temperature, expected_temperature)
    assert f"{shell_count} shells ({shell_count} variables" in caplog.text
    assert "fixed T_e/T_rad=0.900" in caplog.text


def test_prescribed_temperature_still_requires_charge_closure() -> None:
    """Reject a stalled charge residual despite a prescribed temperature."""

    class StalledChargeEvaluator(AnalyticEvaluator):
        def evaluate(
            self,
            electron_density: np.ndarray,
            electron_temperature: np.ndarray,
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            evaluation = super().evaluate(
                electron_density, electron_temperature, level_seed
            )
            evaluation.electron_residual[:] = 1.0e-2
            evaluation.fractional_heating = None
            return evaluation

    with pytest.raises(
        RuntimeError, match="did not close its normalized residuals"
    ):
        PlasmaEquilibriumSolver(fixed_temperature_ratio=0.9).solve(
            StalledChargeEvaluator(),
            np.array([1.0e9]),
            np.array([6000.0]),
            pd.DataFrame([[0.8], [0.2]], columns=[17]),
            np.array([1.0e10]),
            np.array([1.0e4]),
        )


@pytest.mark.parametrize("fixed_temperature_ratio", [None, 0.9])
def test_terminal_refinement_reoptimizes_changed_outer_root(
    fixed_temperature_ratio: float | None,
) -> None:
    """Resume the same outer objective after a terminal context changes it.

    This is orchestration verification: refinement moves a prescribed charge
    root from 2e9 to 3e9 cm^-3 and a thermal root from 8000 to 8500 K.
    The accepted arrays must match the exact evaluation producing the terminal
    diagnostics. Kinetic correctness is tested with the real three-stage atom.
    """

    class RefiningEvaluator:
        level_residual_tolerance = 1e-8

        def __init__(self) -> None:
            self.refined = False
            self.refinements = 0

        def evaluate(
            self,
            electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            density_root = 3e9 if self.refined else 2e9
            temperature_root = 8500.0 if self.refined else 8000.0
            return SimpleNamespace(
                electron_residual=pd.Series(
                    electron_density / density_root - 1
                ),
                fractional_heating=pd.Series(
                    electron_temperature / temperature_root - 1
                ),
                normalized_population=level_seed,
                relative_level_residual=pd.DataFrame(
                    [[0.0 if self.refined else 0.1]]
                ),
                density=electron_density.copy(),
                temperature=electron_temperature.copy(),
            )

        def refine_terminal(
            self,
            electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            previous: SimpleNamespace,
        ) -> SimpleNamespace:
            self.refined = True
            self.refinements += 1
            return self.evaluate(
                electron_density,
                electron_temperature,
                previous.normalized_population,
            )

    evaluator = RefiningEvaluator()
    solution = PlasmaEquilibriumSolver(
        fixed_temperature_ratio=fixed_temperature_ratio
    ).solve(
        evaluator,
        [1e9],
        [7000.0],
        pd.DataFrame([[1.0]]),
        [1e10],
        [1e4],
    )
    npt.assert_allclose(solution.electron_density, [3e9], rtol=1e-10)
    npt.assert_allclose(
        solution.electron_temperature,
        [8500.0 if fixed_temperature_ratio is None else 9000.0],
        rtol=1e-10,
    )
    npt.assert_array_equal(
        solution.electron_density, solution.evaluation.density
    )
    npt.assert_array_equal(
        solution.electron_temperature, solution.evaluation.temperature
    )
    assert evaluator.refinements == 1


def test_terminal_refinement_rejects_unrepairable_levels_after_finite_passes() -> (
    None
):
    """Three ineffective context refinements cannot bypass terminal closure."""

    class UnrepairableEvaluator:
        level_residual_tolerance = 1e-8

        def __init__(self) -> None:
            self.refinements = 0

        def evaluate(
            self,
            electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            del electron_temperature
            return SimpleNamespace(
                electron_residual=pd.Series(np.zeros_like(electron_density)),
                normalized_population=level_seed,
                relative_level_residual=pd.DataFrame([[0.2]]),
            )

        def refine_terminal(
            self,
            electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            previous: SimpleNamespace,
        ) -> SimpleNamespace:
            self.refinements += 1
            return self.evaluate(
                electron_density,
                electron_temperature,
                previous.normalized_population,
            )

    evaluator = UnrepairableEvaluator()
    with pytest.raises(RuntimeError, match="relative level residual"):
        PlasmaEquilibriumSolver(fixed_temperature_ratio=0.9).solve(
            evaluator,
            [1e9],
            [9000.0],
            pd.DataFrame([[1.0]]),
            [1e10],
            [1e4],
        )
    assert evaluator.refinements == 3


def test_terminal_validation_checks_reconstructed_charge() -> None:
    """A density root cannot hide non-neutral returned ion populations.

    This is acceptance wiring: both outer equations vanish, but the supplied
    normalized nuclear charge discrepancy is 1e-3, exceeding 1e-7.
    """

    class NonNeutralEvaluator(AnalyticEvaluator):
        def evaluate(
            self,
            electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            result = super().evaluate(
                electron_density, electron_temperature, level_seed
            )
            result.charge_residual = pd.Series([1e-3])
            return result

    with pytest.raises(
        RuntimeError, match="did not close its normalized residuals"
    ):
        PlasmaEquilibriumSolver().solve(
            NonNeutralEvaluator(),
            [1e9],
            [7000.0],
            pd.DataFrame([[1.0]]),
            [1e10],
            [1e4],
        )


def test_equilibrium_solver_reports_residual_evaluation_progress(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Report nonlinear residual evaluations without changing the solution.

    Claim
    -----
    A long equilibrium solve retains per-evaluation detail at DEBUG while
    keeping INFO for accepted optimizer progress.

    Regime
    ------
    One shell has independent linear density and temperature residuals.

    Verification
    ------------
    The test checks orchestration messages only; the separate analytic-root
    test verifies the numerical result.
    """
    caplog.set_level("DEBUG", logger="tardis.plasma.equilibrium.solver")
    evaluator = AnalyticEvaluator()

    PlasmaEquilibriumSolver().solve(
        evaluator,
        np.array([1.0e9]),
        np.array([6000.0]),
        pd.DataFrame([[0.8], [0.2]], columns=[0]),
        np.array([1.0e10]),
        np.array([1.0e4]),
    )

    assert "Starting plasma equilibrium residual evaluation 1" in caplog.text
    assert "Finished plasma equilibrium residual evaluation 1" in caplog.text
    assert not any(
        record.levelname == "INFO"
        and "plasma equilibrium residual evaluation" in record.message
        for record in caplog.records
    )
    assert (
        "Finished plasma equilibrium least-squares optimization" in caplog.text
    )
    assert "Finished plasma equilibrium terminal validation" in caplog.text


def test_equilibrium_solver_reports_dominant_shells_and_optimizer_cost(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Identify the separate charge and thermal blockers and accepted steps.

    Claim: INFO diagnostics distinguish the two shell-local residuals from
    least-squares progress. Regime: two independent linear shells with a
    known root. Verification: the initial residuals and cost follow directly
    from the evaluator's analytic formulas.
    """
    caplog.set_level("DEBUG", logger="tardis.plasma.equilibrium.solver")
    shell_columns = pd.Index([11, 17])
    level_seed = pd.DataFrame([[0.8, 0.8], [0.2, 0.2]], columns=shell_columns)

    PlasmaEquilibriumSolver().solve(
        AnalyticEvaluator(),
        np.array([1.0e9, 4.0e9]),
        np.array([12000.0, 6000.0]),
        level_seed,
        np.full(2, 1.0e10),
        np.full(2, 1.0e4),
    )

    initial_log = next(
        record.message
        for record in caplog.records
        if "Finished plasma equilibrium residual evaluation 1 "
        in record.message
    )
    assert "electron-density max 1.000000e+00 in shell 17" in initial_log
    assert "trial n_e 4.000000e+09 cm^-3" in initial_log
    assert "thermal max 5.000000e-01 in shell 11" in initial_log
    assert "trial T_e 1.200000e+04 K" in initial_log
    assert "trial cost 7.812500e-01" in initial_log
    first_step_log = next(
        record.message
        for record in caplog.records
        if "Outer optimizer iteration 1" in record.message
    )
    progress = re.search(
        r"cost reduction ([0-9.e+-]+), step infinity norm ([0-9.e+-]+)",
        first_step_log,
    )
    assert progress is not None
    assert float(progress.group(1)) > 0.0
    assert float(progress.group(2)) > 0.0


@pytest.mark.parametrize("heating_offset", [0.0, 2.0])
def test_equilibrium_solver_reports_accepted_thermal_budget(
    caplog: pytest.LogCaptureFixture,
    heating_offset: float,
) -> None:
    """Report the accepted state rather than the last Jacobian perturbation.

    Claim: the displayed heating budget belongs to the accepted thermal
    state. Regime: one shell with a linear thermal root and a temperature-
    dependent synthetic heating budget. Verification: compare the displayed
    heating with the accepted normalized temperature in the same log record.
    """

    class BudgetEvaluator(AnalyticEvaluator):
        def evaluate(
            self,
            electron_density: np.ndarray,
            electron_temperature: np.ndarray,
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            evaluation = super().evaluate(
                electron_density, electron_temperature, level_seed
            )
            evaluation.thermal_budget = pd.DataFrame(
                {
                    "total_heating": electron_temperature / 1000.0
                    + heating_offset,
                    "total_cooling": np.full(len(electron_temperature), 8.0),
                    "bound_free_heating": np.zeros(len(electron_temperature)),
                    "free_bound_cooling": np.zeros(len(electron_temperature)),
                    "free_free_heating": np.zeros(len(electron_temperature)),
                    "free_free_cooling": np.zeros(len(electron_temperature)),
                    "collisional_ionization_heating": np.zeros(
                        len(electron_temperature)
                    ),
                    "collisional_ionization_cooling": np.zeros(
                        len(electron_temperature)
                    ),
                    "collisional_bound_heating": np.zeros(
                        len(electron_temperature)
                    ),
                    "collisional_bound_cooling": np.zeros(
                        len(electron_temperature)
                    ),
                    "nonthermal_heating": np.zeros(len(electron_temperature)),
                    "electron_density": electron_density,
                    "electron_temperature": electron_temperature,
                    "time_explosion": np.full(
                        len(electron_temperature), 172800.0
                    ),
                    "adiabatic_cooling_estimate": np.full(
                        len(electron_temperature), 2.0
                    ),
                },
                index=level_seed.columns,
            )
            return evaluation

    caplog.set_level("INFO", logger="tardis.plasma.equilibrium.solver")
    PlasmaEquilibriumSolver().solve(
        BudgetEvaluator(),
        np.array([1.0e9]),
        np.array([6000.0]),
        pd.DataFrame([[0.8], [0.2]], columns=[17]),
        np.array([1.0e10]),
        np.array([1.0e4]),
    )
    step_logs = [
        record.message
        for record in caplog.records
        if "Outer optimizer iteration" in record.message
    ]
    assert step_logs
    for message in step_logs:
        ratio = float(re.search(r"T_e/T_rad=([0-9.e+-]+)", message).group(1))
        heating = float(re.search(r"H=([0-9.e+-]+)", message).group(1))
        assert heating == pytest.approx(10.0 * ratio + heating_offset, rel=1e-3)
        assert (
            "electron adiabatic estimate (excluded from C)=2.000e+00" in message
        )
        temperature = float(re.search(r"T_e=([0-9.e+-]+) K", message).group(1))
        assert temperature == pytest.approx(1.0e4 * ratio, rel=1e-3)
        assert "t_exp=1.728e+05 s" in message
        deficit = float(re.search(r"H-C=([0-9.e+-]+)", message).group(1))
        assert deficit == pytest.approx(heating - 8.0, abs=0.01)
        if heating_offset == 0.0:
            assert "adiabatic/(H-C)=n/a (H<=C)" in message
        else:
            cooling_ratio = float(
                re.search(
                    r"adiabatic/\(H-C\)=([0-9.]+e[+-][0-9]+)", message
                ).group(1)
            )
            assert cooling_ratio == pytest.approx(2.0 / deficit, rel=1e-3)


def test_equilibrium_solver_supplies_shell_local_jacobian_sparsity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Restrict finite differences to independent shell-local variables.

    Claim
    -----
    Each shell's charge and thermal residuals depend only on that shell's
    electron density and electron temperature.

    Regime
    ------
    Three shells ordered as ``(density, temperature)`` variable pairs.

    Verification
    ------------
    The exact sparsity matrix passed to SciPy contains three independent
    two-by-two diagonal blocks and no cross-shell entries.
    """
    captured_arguments = {}

    def capture_least_squares(
        residual: object,
        initial_guess: np.ndarray,
        **kwargs: object,
    ) -> SimpleNamespace:
        del residual
        captured_arguments.update(kwargs)
        return SimpleNamespace(
            success=True,
            message="test solution",
            x=initial_guess,
            nfev=1,
            njev=1,
            cost=0.0,
            status=1,
        )

    monkeypatch.setattr(
        "tardis.plasma.equilibrium.solver.least_squares",
        capture_least_squares,
    )
    shell_columns = pd.RangeIndex(3)
    level_seed = pd.DataFrame(np.tile([[0.8], [0.2]], 3), columns=shell_columns)

    PlasmaEquilibriumSolver().solve(
        AnalyticEvaluator(),
        np.full(3, 2.0e9),
        np.full(3, 8000.0),
        level_seed,
        np.full(3, 1.0e10),
        np.full(3, 1.0e4),
    )

    expected_sparsity = np.zeros((6, 6), dtype=bool)
    for shell_idx in range(3):
        variable_slice = slice(2 * shell_idx, 2 * shell_idx + 2)
        expected_sparsity[variable_slice, variable_slice] = True
    npt.assert_array_equal(
        captured_arguments["jac_sparsity"].toarray(), expected_sparsity
    )
