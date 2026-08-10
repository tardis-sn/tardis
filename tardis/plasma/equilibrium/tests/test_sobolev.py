import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.opacities.sobolev import calculate_beta_sobolev
from tardis.plasma.equilibrium.population_state import PopulationState
from tardis.plasma.equilibrium.sobolev_solver import SobolevPopulationSolver
from tardis.plasma.exceptions import PlasmaIonizationError

LEVEL_INDEX = pd.MultiIndex.from_tuples(
    [(1, 0, 0), (1, 0, 1)],
    names=["atomic_number", "ion_number", "level_number"],
)
LINE_INDEX = pd.MultiIndex.from_tuples(
    [(1, 0, 0, 1)],
    names=[
        "atomic_number",
        "ion_number",
        "level_number_lower",
        "level_number_upper",
    ],
)
LINES = pd.DataFrame(
    {"wavelength_cm": [1e-5], "f_lu": [1000.0]}, index=LINE_INDEX
)
G = pd.Series([1.0, 1.0], index=LEVEL_INDEX)
METASTABILITY = pd.Series([False, False], index=LEVEL_INDEX)


class BetaPopulationSolver:
    """Return a deterministic population state whose excitation tracks beta."""

    @staticmethod
    def _state(beta_sobolev: pd.DataFrame) -> PopulationState:
        """Build the deterministic state for one beta value."""
        beta_value = beta_sobolev.iloc[0, 0]
        level_number_density = pd.DataFrame(
            [[1.0], [0.2 + 0.4 * beta_value]],
            index=LEVEL_INDEX,
            columns=[0],
        )
        ion_number_density = pd.DataFrame(
            [[1.0 + beta_value]],
            index=pd.MultiIndex.from_tuples(
                [(1, 0)], names=["atomic_number", "ion_number"]
            ),
            columns=[0],
        )
        return PopulationState(
            electron_densities=pd.Series([beta_value], index=[0]),
            elemental_populations={},
            ion_number_density=ion_number_density,
            level_number_density=level_number_density,
            level_boltzmann_factor=level_number_density.copy(),
        )

    def solve(
        self,
        beta_sobolev: pd.DataFrame,
        electron_density_seed: pd.Series | None = None,
    ) -> PopulationState:
        del electron_density_seed
        return self._state(beta_sobolev)


def test_sobolev_fixed_point_converges_from_positive_states() -> None:
    """Distinct positive seeds must converge to one population/Sobolev state."""
    results = []
    for upper_population in (0.1, 0.8):
        initial_level_number_density = pd.DataFrame(
            [[1.0], [upper_population]], index=LEVEL_INDEX, columns=[0]
        )
        solver = SobolevPopulationSolver(
            BetaPopulationSolver(),
            LINES,
            1 * u.day,
            G,
            [0],
            [1],
            METASTABILITY,
            tolerance=1e-7,
        )
        results.append(solver.solve(initial_level_number_density))

    first_state, _, first_tau, first_beta = results[0]
    second_state, _, _, second_beta = results[1]
    assert np.all(first_state.level_number_density.to_numpy() > 0.0)
    pd.testing.assert_frame_equal(
        first_state.level_number_density,
        second_state.level_number_density,
        rtol=1e-7,
        atol=0.0,
    )
    pd.testing.assert_frame_equal(
        first_beta,
        second_beta,
        rtol=1e-7,
        atol=0.0,
    )
    np.testing.assert_allclose(
        first_state.electron_densities, first_beta.iloc[0]
    )
    np.testing.assert_allclose(
        first_beta.to_numpy(),
        calculate_beta_sobolev(first_tau).to_numpy(),
    )


def test_sobolev_nonconvergence_reports_final_update_norms() -> None:
    """A failed population/escape-probability closure must expose its residuals."""
    initial_level_number_density = pd.DataFrame(
        [[1.0], [0.1]], index=LEVEL_INDEX, columns=[0]
    )
    solver = SobolevPopulationSolver(
        BetaPopulationSolver(),
        LINES,
        1 * u.day,
        G,
        [0],
        [1],
        METASTABILITY,
        max_iterations=1,
    )

    with pytest.raises(
        PlasmaIonizationError,
        match=r"update_norm=.*population_norm=.*relaxation=",
    ):
        solver.solve(initial_level_number_density)


def test_sobolev_populations_respond_to_optical_depth_scale() -> None:
    """Explosion time must change both converged escape and excitation states."""
    initial_level_number_density = pd.DataFrame(
        [[1.0], [0.1]], index=LEVEL_INDEX, columns=[0]
    )
    results = []
    for time_explosion in (1.0, 2.0):
        solver = SobolevPopulationSolver(
            BetaPopulationSolver(),
            LINES,
            time_explosion * u.day,
            G,
            [0],
            [1],
            METASTABILITY,
            tolerance=1e-7,
        )
        results.append(solver.solve(initial_level_number_density))

    first_state, _, _, first_beta = results[0]
    second_state, _, _, second_beta = results[1]
    assert not np.allclose(first_beta, second_beta)
    assert not np.allclose(
        first_state.level_number_density,
        second_state.level_number_density,
    )


def test_sobolev_retries_a_growing_update_with_half_relaxation() -> None:
    """A growing fixed-point residual must retry the previous update halfway."""

    class RetryPopulationSolver:
        """Return populations that make the first undamped update overshoot."""

        def __init__(self) -> None:
            self.beta_values: list[float] = []

        def solve(
            self,
            beta_sobolev: pd.DataFrame,
            electron_density_seed: pd.Series | None = None,
        ) -> PopulationState:
            del electron_density_seed
            beta_value = float(beta_sobolev.iloc[0, 0])
            self.beta_values.append(beta_value)
            if np.isclose(beta_value, 0.2):
                upper_population = 0.8
            elif np.isclose(beta_value, 0.8):
                upper_population = 0.0
            else:
                upper_population = 0.5
            level_number_density = pd.DataFrame(
                [[1.0], [upper_population]],
                index=LEVEL_INDEX,
                columns=[0],
            )
            return PopulationState(
                electron_densities=pd.Series([1.0], index=[0]),
                elemental_populations={},
                ion_number_density=pd.DataFrame(
                    [[1.0]],
                    index=pd.MultiIndex.from_tuples(
                        [(1, 0)], names=["atomic_number", "ion_number"]
                    ),
                    columns=[0],
                ),
                level_number_density=level_number_density,
                level_boltzmann_factor=level_number_density.copy(),
            )

    class DirectBetaSobolevSolver(SobolevPopulationSolver):
        """Use the upper-level population directly as beta for this test."""

        def _calculate_sobolev_state(
            self, level_number_density: pd.DataFrame
        ) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
            beta_sobolev = pd.DataFrame(
                [[level_number_density.iloc[1, 0]]],
                index=LINE_INDEX,
                columns=[0],
            )
            return np.ones((1, 1)), beta_sobolev, beta_sobolev

    population_solver = RetryPopulationSolver()
    solver = DirectBetaSobolevSolver(
        population_solver,
        LINES,
        1 * u.day,
        G,
        [0],
        [1],
        METASTABILITY,
        tolerance=1e-7,
    )
    initial_level_number_density = pd.DataFrame(
        [[1.0], [0.2]], index=LEVEL_INDEX, columns=[0]
    )
    initial_beta_sobolev = pd.DataFrame([[0.2]], index=LINE_INDEX, columns=[0])

    state, _, _, beta_sobolev = solver.solve(
        initial_level_number_density,
        initial_beta_sobolev,
    )

    np.testing.assert_allclose(
        population_solver.beta_values[:3], [0.2, 0.8, 0.5]
    )
    np.testing.assert_allclose(state.level_number_density.iloc[1], 0.5)
    np.testing.assert_allclose(beta_sobolev, 0.5)


def test_sobolev_rejects_final_population_drift() -> None:
    """Final back-substitution must preserve populations as well as beta."""
    initial_level_number_density = pd.DataFrame(
        [[1.0], [0.1]], index=LEVEL_INDEX, columns=[0]
    )
    baseline_solver = SobolevPopulationSolver(
        BetaPopulationSolver(),
        LINES,
        1 * u.day,
        G,
        [0],
        [1],
        METASTABILITY,
        tolerance=1e-7,
    )
    baseline_state, _, _, baseline_beta = baseline_solver.solve(
        initial_level_number_density
    )

    class FinalDriftPopulationSolver(BetaPopulationSolver):
        """Change populations only during the final back-substitution."""

        def __init__(self) -> None:
            self.calls = 0

        def solve(
            self,
            beta_sobolev: pd.DataFrame,
            electron_density_seed: pd.Series | None = None,
        ) -> PopulationState:
            state = super().solve(beta_sobolev, electron_density_seed)
            self.calls += 1
            if self.calls != 2:
                return state
            level_number_density = state.level_number_density + 1.0
            return PopulationState(
                electron_densities=state.electron_densities,
                elemental_populations=state.elemental_populations,
                ion_number_density=state.ion_number_density,
                level_number_density=level_number_density,
                level_boltzmann_factor=level_number_density.copy(),
            )

    solver = SobolevPopulationSolver(
        FinalDriftPopulationSolver(),
        LINES,
        1 * u.day,
        G,
        [0],
        [1],
        METASTABILITY,
        tolerance=1e-7,
    )

    with pytest.raises(PlasmaIonizationError, match="population_norm"):
        solver.solve(
            baseline_state.level_number_density,
            baseline_beta,
        )
