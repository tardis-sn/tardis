import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest
from astropy import units as u

from tardis.plasma.equilibrium.rates.collision_strengths import (
    UpsilonRegemorterSolver,
    exp1_times_exp,
)

EXPECTED_INDEX_NAMES = [
    "atomic_number",
    "ion_number",
    "level_number_lower",
    "level_number_upper",
]


def _make_transition_data(
    atomic_number=1,
    ion_number=0,
    level_lower=0,
    level_upper=1,
    f_lu=0.4162,
    nu=2.47e15,
):
    """Build a minimal valid transition DataFrame with one transition."""
    index = pd.MultiIndex.from_tuples(
        [(atomic_number, ion_number, level_lower, level_upper)],
        names=EXPECTED_INDEX_NAMES,
    )
    return pd.DataFrame({"f_lu": [f_lu], "nu": [nu]}, index=index)


@pytest.fixture
def minimal_transition_data():
    """Single-transition DataFrame for Lyman-alpha (H I 1->2)."""
    return _make_transition_data()


@pytest.fixture
def minimal_solver(minimal_transition_data):
    """UpsilonRegemorterSolver built from minimal single-transition data."""
    return UpsilonRegemorterSolver(minimal_transition_data)


@pytest.fixture
def sample_temperatures():
    """Representative electron temperatures spanning typical nebular conditions."""
    return np.array([5000.0, 10000.0, 20000.0]) * u.K


class TestExp1TimesExp:
    """Tests for the exp1_times_exp utility function."""

    def test_returns_array_of_same_shape(self):
        """Output shape matches input shape."""
        x = np.array([1.0, 2.0, 3.0])
        result = exp1_times_exp(x)
        assert result.shape == x.shape

    def test_small_values_close_to_scipy(self):
        """For small x, exp1(x)*exp(x) matches the direct scipy computation."""
        from scipy.special import exp1

        x = np.array([0.1, 0.5, 1.0, 5.0])
        expected = exp1(x) * np.exp(x)
        result = exp1_times_exp(x)
        npt.assert_allclose(result, expected, rtol=1e-12)

    def test_large_values_use_laurent_series(self):
        """For x > 500, function uses Laurent series; result must be finite and positive."""
        x = np.array([501.0, 1000.0, 5000.0])
        result = exp1_times_exp(x)
        assert np.all(np.isfinite(result))
        assert np.all(result > 0)

    def test_large_values_match_asymptotic_approximation(self):
        """Laurent series result is close to 1/x for very large x."""
        x = np.array([1000.0, 10000.0])
        result = exp1_times_exp(x)
        # leading term of Laurent series is x^-1
        npt.assert_allclose(result, 1.0 / x, rtol=1e-2)

    def test_monotonically_decreasing(self):
        """exp1(x)*exp(x) is monotonically decreasing for positive x."""
        x = np.linspace(0.1, 10.0, 50)
        result = exp1_times_exp(x)
        assert np.all(np.diff(result) < 0)


class TestUpsilonRegemorterSolverInit:
    """Tests for __init__ validation logic."""

    def test_init_valid(self, minimal_solver):
        """Solver constructs without error from valid transition data."""
        assert minimal_solver is not None
        assert hasattr(minimal_solver, "transition_data")

    def test_init_default_g_bar(self, minimal_solver):
        """Default g_bar is 0.2 as per van Regemorter approximation."""
        assert minimal_solver.g_bar == 0.2

    def test_init_custom_g_bar(self, minimal_transition_data):
        """Custom g_bar is stored correctly."""
        solver = UpsilonRegemorterSolver(minimal_transition_data, g_bar=0.7)
        assert solver.g_bar == 0.7

    def test_init_missing_f_lu_raises(self, minimal_transition_data):
        """AssertionError raised when f_lu column is missing."""
        df = minimal_transition_data.drop(columns=["f_lu"])
        with pytest.raises(AssertionError):
            UpsilonRegemorterSolver(df)

    def test_init_missing_nu_raises(self, minimal_transition_data):
        """AssertionError raised when nu column is missing."""
        df = minimal_transition_data.drop(columns=["nu"])
        with pytest.raises(AssertionError):
            UpsilonRegemorterSolver(df)

    def test_init_wrong_index_names_raises(self, minimal_transition_data):
        """AssertionError raised when index level names are wrong."""
        df = minimal_transition_data.copy()
        df.index.names = ["atomic_number", "ion_number", "lower", "upper"]
        with pytest.raises(AssertionError):
            UpsilonRegemorterSolver(df)

    def test_init_inverted_level_numbers_raises(self):
        """AssertionError raised when level_number_lower >= level_number_upper."""
        index = pd.MultiIndex.from_tuples(
            [(1, 0, 1, 0)], names=EXPECTED_INDEX_NAMES
        )
        df = pd.DataFrame({"f_lu": [0.4], "nu": [2.47e15]}, index=index)
        with pytest.raises(AssertionError):
            UpsilonRegemorterSolver(df)

    def test_init_equal_level_numbers_raises(self):
        """AssertionError raised when level_number_lower == level_number_upper."""
        index = pd.MultiIndex.from_tuples(
            [(1, 0, 0, 0)], names=EXPECTED_INDEX_NAMES
        )
        df = pd.DataFrame({"f_lu": [0.4], "nu": [2.47e15]}, index=index)
        with pytest.raises(AssertionError):
            UpsilonRegemorterSolver(df)

    def test_init_sorts_index(self, minimal_transition_data):
        """Transition data is sorted by index after construction."""
        extra_index = pd.MultiIndex.from_tuples(
            [(1, 0, 2, 3)], names=EXPECTED_INDEX_NAMES
        )
        extra_row = pd.DataFrame(
            {"f_lu": [0.1], "nu": [5e14]}, index=extra_index
        )
        unsorted = pd.concat([extra_row, minimal_transition_data])
        solver = UpsilonRegemorterSolver(unsorted)
        assert solver.transition_data.index.is_monotonic_increasing


class TestUpsilonRegemorterSolverOutputStructure:
    """Tests for the structure of the DataFrame returned by solve()."""

    def test_solve_returns_dataframe(self, minimal_solver, sample_temperatures):
        """solve() returns a pandas DataFrame."""
        result = minimal_solver.solve(sample_temperatures)
        assert isinstance(result, pd.DataFrame)

    def test_solve_output_index_names(self, minimal_solver, sample_temperatures):
        """Output index has the same 4-level MultiIndex as the input."""
        result = minimal_solver.solve(sample_temperatures)
        assert list(result.index.names) == EXPECTED_INDEX_NAMES

    def test_solve_output_row_count(self, minimal_solver, sample_temperatures):
        """Number of output rows equals number of input transitions."""
        result = minimal_solver.solve(sample_temperatures)
        assert len(result) == len(minimal_solver.transition_data)

    def test_solve_output_column_count(self, minimal_solver, sample_temperatures):
        """Number of output columns equals number of temperature points."""
        result = minimal_solver.solve(sample_temperatures)
        assert result.shape[1] == len(sample_temperatures)

    def test_solve_multi_transition(self, sample_temperatures):
        """solve() handles multiple transitions without error."""
        idx = pd.MultiIndex.from_tuples(
            [(1, 0, 0, 1), (1, 0, 0, 2), (1, 0, 1, 2)],
            names=EXPECTED_INDEX_NAMES,
        )
        df = pd.DataFrame(
            {"f_lu": [0.42, 0.08, 0.01], "nu": [2.47e15, 1.6e15, 8e14]},
            index=idx,
        )
        solver = UpsilonRegemorterSolver(df)
        result = solver.solve(sample_temperatures)
        assert result.shape == (3, 3)


class TestUpsilonRegemorterSolverPhysics:
    """Tests verifying physical properties of the van Regemorter collision strengths."""

    def test_collision_strengths_positive(self, minimal_solver, sample_temperatures):
        """All upsilon_g values must be positive for physical inputs."""
        result = minimal_solver.solve(sample_temperatures)
        assert np.all(result.values > 0)

    def test_strength_temperature_independent_in_gbar_regime(
        self, minimal_transition_data
    ):
        """Upsilon_g is temperature-independent when g_bar dominates.

        In the van Regemorter formula, T * u0 = h*nu / k_B is a constant,
        so when gamma = g_bar (constant), the result does not depend on T.
        """
        low_T = np.array([5000.0]) * u.K
        high_T = np.array([20000.0]) * u.K
        solver = UpsilonRegemorterSolver(minimal_transition_data)

        result_low = solver.solve(low_T)
        result_high = solver.solve(high_T)

        npt.assert_allclose(
            result_high.values[0, 0], result_low.values[0, 0], rtol=1e-10
        )

    def test_larger_f_lu_gives_larger_strength(self, sample_temperatures):
        """Higher oscillator strength f_lu leads to larger collision strength."""
        df_low = _make_transition_data(f_lu=0.1)
        df_high = _make_transition_data(f_lu=0.8)

        result_low = UpsilonRegemorterSolver(df_low).solve(sample_temperatures)
        result_high = UpsilonRegemorterSolver(df_high).solve(sample_temperatures)

        assert np.all(result_high.values > result_low.values)

    def test_g_bar_07_gives_larger_strength_than_02(self, minimal_transition_data, sample_temperatures):
        """g_bar=0.7 yields larger collision strengths than g_bar=0.2."""
        solver_low_g = UpsilonRegemorterSolver(minimal_transition_data, g_bar=0.2)
        solver_high_g = UpsilonRegemorterSolver(minimal_transition_data, g_bar=0.7)

        result_low = solver_low_g.solve(sample_temperatures)
        result_high = solver_high_g.solve(sample_temperatures)

        # At low u0 (high T), gamma = max(g_bar, component) may be dominated by
        # the gamma_component, so check only the low-temperature shell where
        # the g_bar term is most likely to dominate.
        assert result_high.values[0, 0] >= result_low.values[0, 0]

    def test_strength_scales_with_oscillator_strength(self, sample_temperatures):
        """Upsilon_g scales linearly with f_lu (van Regemorter approximation)."""
        df_1 = _make_transition_data(f_lu=0.2)
        df_2 = _make_transition_data(f_lu=0.4)

        result_1 = UpsilonRegemorterSolver(df_1).solve(sample_temperatures)
        result_2 = UpsilonRegemorterSolver(df_2).solve(sample_temperatures)

        npt.assert_allclose(result_2.values, 2 * result_1.values, rtol=1e-10)
