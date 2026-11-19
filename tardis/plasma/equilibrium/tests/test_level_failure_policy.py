"""Verify narrow continuation of finite terminal kinetic-closure failure."""

from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest

from tardis.plasma.equilibrium.solver import PlasmaEquilibriumSolver


@pytest.fixture
def terminal_policy_evaluator() -> object:
    """Supply an analytic charge root and conservative two-level helium state."""

    class Evaluator:
        level_residual_tolerance = 1e-8
        equilibrium_species = (2, 0)
        elemental_number_density = pd.DataFrame(
            [[1e9]], index=[2], columns=[17]
        )

        def __init__(self) -> None:
            self.refinements = 0
            self.level_residual = -0.1921
            self.invalid_field: tuple[str, float] | None = None
            self.electron_offset = 0.0
            self.thermal_offset = 0.0
            self.charge_offset = 0.0
            self.corruption: str | None = None
            self.close_after: int | None = None
            self.invalid_density: tuple[str, float] | None = None
            self.subnormal_parent: tuple[int, bool] | None = None

        def evaluate(
            self,
            electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            level_seed: pd.DataFrame,
        ) -> SimpleNamespace:
            ion_index = pd.MultiIndex.from_tuples(
                [(2, 0), (2, 1)], names=["atomic_number", "ion_number"]
            )
            level_index = pd.MultiIndex.from_tuples(
                [(2, 0, 3), (2, 0, 9)],
                names=["atomic_number", "ion_number", "level_number"],
            )
            result = SimpleNamespace(
                electron_residual=pd.Series(
                    electron_density / 4e8 - 1 + self.electron_offset,
                    index=[17],
                ),
                fractional_heating=pd.Series(
                    electron_temperature / 8000 - 1 + self.thermal_offset,
                    index=[17],
                ),
                charge_residual=pd.Series([self.charge_offset], index=[17]),
                normalized_population=level_seed.copy().set_axis(level_index),
                ion_population=pd.DataFrame(
                    [[6e8], [4e8]], index=ion_index, columns=[17]
                ),
                absolute_level_population=pd.DataFrame(
                    [[4.5e8], [1.5e8]], index=level_index, columns=[17]
                ),
                tau_sobolev=pd.DataFrame([[-1.0]], columns=[17]),
                beta_sobolev=pd.DataFrame([[1.5]], columns=[17]),
                relative_level_residual=pd.DataFrame(
                    [[self.level_residual], [0.0]],
                    index=level_index,
                    columns=[17],
                ),
                density=electron_density.copy(),
                temperature=electron_temperature.copy(),
            )
            if self.invalid_density is not None:
                field, value = self.invalid_density
                setattr(result, field, pd.Series([value], index=[17]))
            if self.subnormal_parent is not None:
                quanta, represented = self.subnormal_parent
                quantum = np.nextafter(0.0, 1.0)
                result.ion_population.loc[(28, 8), :] = quanta * quantum
                for level_number in range(26):
                    result.absolute_level_population.loc[
                        (28, 8, level_number), :
                    ] = quantum if represented else 0.0
            if self.invalid_field is not None:
                field, value = self.invalid_field
                getattr(result, field).iloc[0, 0] = value
            if self.corruption == "fractions":
                result.normalized_population.iloc[:, 0] = [0.25, 0.75]
            elif self.corruption == "shells":
                result.ion_population.columns = [18]
            elif self.corruption == "missing_parent":
                result.ion_population = result.ion_population.drop((2, 0))
            elif self.corruption == "missing_diagnostics":
                del result.relative_level_residual
            return result

        def refine_terminal(
            self,
            electron_density: npt.NDArray[np.float64],
            electron_temperature: npt.NDArray[np.float64],
            previous: SimpleNamespace,
        ) -> SimpleNamespace:
            self.refinements += 1
            if self.close_after == self.refinements:
                self.level_residual = 0.0
            return self.evaluate(
                electron_density,
                electron_temperature,
                previous.normalized_population,
            )

    return Evaluator()


@pytest.mark.parametrize("policy", ["raise", "warn"])
def test_terminal_level_policy_retains_exact_candidate(
    terminal_policy_evaluator: object,
    policy: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Finite open kinetics may continue, never claim physical equilibrium.

    Regime: analytic n_e=4e8 charge root, prescribed temperature, normalized
    [.75,.25] levels and conserved 1e9 helium nuclei. The explicit signed
    residual -.1921 remains above the unchanged 1e-8 target after three
    refinements. This is acceptance orchestration, not kinetic validation.
    """
    solver = PlasmaEquilibriumSolver(
        fixed_temperature_ratio=0.9, level_failure_policy=policy
    )
    arguments = (
        terminal_policy_evaluator,
        [2e8],
        [7000.0],
        pd.DataFrame([[0.75], [0.25]], columns=[17]),
        [2e9],
        [1e4],
    )
    if policy == "raise":
        with pytest.raises(RuntimeError, match="relative level residual"):
            solver.solve(*arguments)
    else:
        solution = solver.solve(*arguments)
        assert solution.level_converged is False
        assert solution.maximum_relative_level_residual == pytest.approx(0.1921)
        npt.assert_allclose(solution.electron_density, [4e8], rtol=1e-10)
        npt.assert_array_equal(
            solution.electron_density, solution.evaluation.density
        )
        npt.assert_array_equal(
            solution.electron_temperature, solution.evaluation.temperature
        )
        assert "NOT converged level equilibrium" in caplog.text
        assert "(2, 0, 3)" in caplog.text and "shell 17" in caplog.text
        assert "-1.921e-01" in caplog.text and "refinements=3" in caplog.text
    assert terminal_policy_evaluator.refinements == 3


@pytest.mark.parametrize("policy", ["raise", "warn"])
@pytest.mark.parametrize("residual", [0.0, 0.5e-8])
def test_closed_level_policy_reports_measured_closure(
    terminal_policy_evaluator: object,
    policy: str,
    residual: float,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Both policies accept the same unchanged kinetic tolerance boundary."""
    terminal_policy_evaluator.level_residual = residual
    solution = PlasmaEquilibriumSolver(
        fixed_temperature_ratio=0.9, level_failure_policy=policy
    ).solve(
        terminal_policy_evaluator,
        [4e8],
        [9000.0],
        pd.DataFrame([[0.75], [0.25]], columns=[17]),
        [2e9],
        [1e4],
    )
    assert solution.level_converged is True
    assert solution.maximum_relative_level_residual == residual
    assert terminal_policy_evaluator.refinements == 0
    assert not [
        record for record in caplog.records if record.levelname == "WARNING"
    ]


@pytest.mark.parametrize(
    "field,value,reason",
    [
        ("relative_level_residual", np.nan, "nonfinite relative level"),
        ("relative_level_residual", np.inf, "nonfinite relative level"),
        ("absolute_level_population", np.nan, "invalid absolute_level"),
        ("absolute_level_population", -1.0, "invalid absolute_level"),
        ("ion_population", -1.0, "invalid ion_population"),
        ("normalized_population", -1.0, "invalid normalized_population"),
        ("normalized_population", 0.5, "normalization"),
        ("absolute_level_population", 5e8, "level/ion totals"),
        ("ion_population", 7e8, "level/ion totals"),
        ("tau_sobolev", np.inf, "invalid tau_sobolev"),
        ("beta_sobolev", np.nan, "invalid beta_sobolev"),
    ],
)
def test_warn_policy_rejects_invalid_terminal_physics(
    terminal_policy_evaluator: object,
    field: str,
    value: float,
    reason: str,
) -> None:
    """Warning is not permission to return invalid or unmatched physical data.

    Explicit violations are orders above binary64 normalization arithmetic;
    no failed physical balance is treated as an acceptable kinetic warning.
    These are boundary-orchestration tests, not helium equilibrium evidence.
    """
    terminal_policy_evaluator.invalid_field = (field, value)
    with pytest.raises(RuntimeError, match=reason):
        PlasmaEquilibriumSolver(
            fixed_temperature_ratio=0.9, level_failure_policy="warn"
        ).solve(
            terminal_policy_evaluator,
            [4e8],
            [9000.0],
            pd.DataFrame([[0.75], [0.25]], columns=[17]),
            [2e9],
            [1e4],
        )


@pytest.mark.parametrize(
    "mode", ["charge", "electron", "thermal", "unused_thermal"]
)
def test_warn_policy_preserves_outer_residual_acceptance(
    terminal_policy_evaluator: object,
    mode: str,
) -> None:
    """Failed reconstructed charge and active thermal closure remain fatal."""
    terminal_policy_evaluator.charge_offset = 1e-3 if mode == "charge" else 0.0
    terminal_policy_evaluator.electron_offset = (
        10.0 if mode == "electron" else 0.0
    )
    terminal_policy_evaluator.thermal_offset = (
        np.nan
        if mode == "unused_thermal"
        else (10.0 if mode == "thermal" else 0.0)
    )
    solver = PlasmaEquilibriumSolver(
        fixed_temperature_ratio=None if mode == "thermal" else 0.9,
        level_failure_policy="warn",
    )
    arguments = (
        terminal_policy_evaluator,
        [4e8],
        [9000.0],
        pd.DataFrame([[0.75], [0.25]], columns=[17]),
        [2e9],
        [1e4],
    )
    if mode == "unused_thermal":
        assert solver.solve(*arguments).level_converged is False
    else:
        with pytest.raises(RuntimeError, match="normalized residuals"):
            solver.solve(*arguments)


def test_level_policy_rejects_unknown_setting() -> None:
    """The explicit user-facing policy setting has an actionable error."""
    with pytest.raises(ValueError, match="level_failure_policy"):
        PlasmaEquilibriumSolver(level_failure_policy="ignore")


@pytest.mark.parametrize(
    "corruption,reason",
    [
        ("fractions", "selected-level fractions"),
        ("shells", "misaligned ion_population"),
        ("missing_parent", "level/ion totals"),
        ("nuclei", "conserved nuclei"),
        ("missing_element", "conserved nuclei"),
        ("extra_element", "conserved nuclei"),
    ],
)
def test_warn_policy_checks_labeled_physical_consistency(
    terminal_policy_evaluator: object,
    corruption: str,
    reason: str,
) -> None:
    """Per-ion totals cannot hide wrong fractions or missing nuclear species."""
    evaluator = terminal_policy_evaluator
    evaluator.corruption = corruption
    evaluator.elemental_number_density = (
        evaluator.elemental_number_density.copy()
    )
    if corruption == "nuclei":
        evaluator.elemental_number_density.loc[2] = 1.1e9
    elif corruption == "missing_element":
        evaluator.elemental_number_density.loc[1] = 1e8
    elif corruption == "extra_element":
        evaluator.elemental_number_density = (
            evaluator.elemental_number_density.drop(2)
        )
    with pytest.raises(RuntimeError, match=reason):
        PlasmaEquilibriumSolver(
            fixed_temperature_ratio=0.9, level_failure_policy="warn"
        ).solve(
            evaluator,
            [4e8],
            [9000.0],
            pd.DataFrame([[0.75], [0.25]], columns=[17]),
            [2e9],
            [1e4],
        )


@pytest.mark.parametrize("policy", ["raise", "warn"])
def test_terminal_level_refinement_recovers_current_status(
    terminal_policy_evaluator: object,
    policy: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A second successful terminal refinement yields true current closure."""
    terminal_policy_evaluator.close_after = 2
    solution = PlasmaEquilibriumSolver(
        fixed_temperature_ratio=0.9, level_failure_policy=policy
    ).solve(
        terminal_policy_evaluator,
        [4e8],
        [9000.0],
        pd.DataFrame([[0.75], [0.25]], columns=[17]),
        [2e9],
        [1e4],
    )
    assert solution.level_converged is True
    assert solution.maximum_relative_level_residual == 0.0
    assert terminal_policy_evaluator.refinements == 2
    assert "NOT converged" not in caplog.text


def test_missing_optional_kinetic_diagnostics_are_not_measured_zero(
    terminal_policy_evaluator: object,
) -> None:
    """Preserve the prior optional evaluator contract without inventing data."""
    terminal_policy_evaluator.corruption = "missing_diagnostics"
    solution = PlasmaEquilibriumSolver(fixed_temperature_ratio=0.9).solve(
        terminal_policy_evaluator,
        [4e8],
        [9000.0],
        pd.DataFrame([[0.75], [0.25]], columns=[17]),
        [2e9],
        [1e4],
    )
    assert solution.level_converged is True
    assert solution.maximum_relative_level_residual is None
    assert terminal_policy_evaluator.refinements == 0


@pytest.mark.parametrize(
    "field", ["trial_electron_density", "charge_solved_electron_density"]
)
@pytest.mark.parametrize("value", [np.nan, np.inf, 0.0, -1.0])
def test_warn_policy_rejects_invalid_optional_physical_density(
    terminal_policy_evaluator: object,
    field: str,
    value: float,
) -> None:
    """A closed mock residual cannot hide an invalid actual density field."""
    terminal_policy_evaluator.invalid_density = (field, value)
    with pytest.raises(RuntimeError, match=f"invalid {field}"):
        PlasmaEquilibriumSolver(
            fixed_temperature_ratio=0.9, level_failure_policy="warn"
        ).solve(
            terminal_policy_evaluator,
            [4e8],
            [9000.0],
            pd.DataFrame([[0.75], [0.25]], columns=[17]),
            [2e9],
            [1e4],
        )


@pytest.mark.parametrize(
    "quanta,represented,accepted",
    [(13, False, True), (14, False, False), (26, True, True)],
)
def test_terminal_physical_guard_preserves_exact_allocation_underflow_contract(
    terminal_policy_evaluator: object,
    quanta: int,
    represented: bool,
    accepted: bool,
) -> None:
    """Only exact all-zero allocation below floor(N/2) quanta may be orphaned.

    Binary64 round-to-nearest allows 26 normalized level shares to all vanish
    for at most 13 smallest-subnormal quanta. One quantum above is not repaired
    or accepted; a represented one-quantum-per-level distribution is retained.
    This verifies representability, not a physical abundance cutoff.
    """
    evaluator = terminal_policy_evaluator
    evaluator.subnormal_parent = (quanta, represented)
    evaluator.elemental_number_density = (
        evaluator.elemental_number_density.copy()
    )
    evaluator.elemental_number_density.loc[28] = quanta * np.nextafter(0.0, 1.0)
    solver = PlasmaEquilibriumSolver(
        fixed_temperature_ratio=0.9, level_failure_policy="warn"
    )
    arguments = (
        evaluator,
        [4e8],
        [9000.0],
        pd.DataFrame([[0.75], [0.25]], columns=[17]),
        [2e9],
        [1e4],
    )
    if accepted:
        result = solver.solve(*arguments)
        npt.assert_equal(
            result.evaluation.ion_population.loc[(28, 8), 17],
            quanta * np.nextafter(0.0, 1.0),
        )
    else:
        with pytest.raises(RuntimeError, match="level/ion totals"):
            solver.solve(*arguments)
