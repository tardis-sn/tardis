import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import pytest
from astropy import units as u

from tardis.plasma.equilibrium.rates import RadiativeRatesSolver
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
    PlanckianRadiationField,
)


# ---------------------------------------------------------------------------
# Helpers / shared mock data
# ---------------------------------------------------------------------------

EXPECTED_INDEX_NAMES = [
    "atomic_number",
    "ion_number",
    "level_number_lower",
    "level_number_upper",
]

EXPECTED_OUTPUT_INDEX_NAMES = [
    "atomic_number",
    "ion_number",
    "ion_number_source",
    "ion_number_destination",
    "level_number_source",
    "level_number_destination",
]


def _make_einstein_coefficients(
    atomic_number=1,
    ion_number=0,
    level_lower=0,
    level_upper=1,
    A_ul=6.265e8,
    nu=2.47e15,
):
    """
    Build a minimal valid Einstein coefficient DataFrame with one transition.

    Uses the Einstein relation:
        B_ul = A_ul * c^2 / (2 * h * nu^3)
        B_lu = g_u / g_l * B_ul   (assuming g_u = g_l = 1 for simplicity)
    """
    from tardis import constants

    c = constants.c.cgs.value
    h = constants.h.cgs.value

    B_ul = A_ul * c**2 / (2 * h * nu**3)
    B_lu = B_ul  # equal statistical weights assumed

    index = pd.MultiIndex.from_tuples(
        [(atomic_number, ion_number, level_lower, level_upper)],
        names=EXPECTED_INDEX_NAMES,
    )
    return pd.DataFrame(
        {"A_ul": [A_ul], "B_ul": [B_ul], "B_lu": [B_lu], "nu": [nu]},
        index=index,
    )


@pytest.fixture
def minimal_einstein_coefficients():
    """Single-transition Einstein coefficient DataFrame."""
    return _make_einstein_coefficients()


@pytest.fixture
def minimal_rad_field():
    """Single-shell DilutePlanckianRadiationField."""
    return DilutePlanckianRadiationField(
        temperature=np.array([10000.0]) * u.K,
        dilution_factor=np.array([0.5]),
    )


@pytest.fixture
def minimal_solver(minimal_einstein_coefficients):
    """RadiativeRatesSolver built from the minimal single-transition data."""
    return RadiativeRatesSolver(minimal_einstein_coefficients)


# ---------------------------------------------------------------------------
# Construction / Validation Tests
# ---------------------------------------------------------------------------


class TestRadiativeRatesSolverInit:
    """Tests for __init__ validation logic."""

    def test_init_valid(self, radiative_rate_solver):
        """Solver constructs without error from real Si Chianti atomic data."""
        assert radiative_rate_solver is not None
        assert hasattr(radiative_rate_solver, "einstein_coefficients")

    def test_init_stores_einstein_coefficients(self, minimal_solver, minimal_einstein_coefficients):
        """Einstein coefficients DataFrame is stored on the object."""
        assert minimal_solver.einstein_coefficients is not None
        assert set(minimal_solver.einstein_coefficients.columns) >= {
            "A_ul", "B_ul", "B_lu", "nu"
        }

    def test_init_missing_A_ul_raises(self, minimal_einstein_coefficients):
        """AssertionError raised when A_ul column is missing."""
        df = minimal_einstein_coefficients.drop(columns=["A_ul"])
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(df)

    def test_init_missing_B_ul_raises(self, minimal_einstein_coefficients):
        """AssertionError raised when B_ul column is missing."""
        df = minimal_einstein_coefficients.drop(columns=["B_ul"])
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(df)

    def test_init_missing_B_lu_raises(self, minimal_einstein_coefficients):
        """AssertionError raised when B_lu column is missing."""
        df = minimal_einstein_coefficients.drop(columns=["B_lu"])
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(df)

    def test_init_missing_nu_raises(self, minimal_einstein_coefficients):
        """AssertionError raised when nu column is missing."""
        df = minimal_einstein_coefficients.drop(columns=["nu"])
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(df)

    def test_init_wrong_index_names_raises(self, minimal_einstein_coefficients):
        """AssertionError raised when index level names are wrong."""
        df = minimal_einstein_coefficients.copy()
        df.index.names = [
            "atomic_number", "ion_number", "level_lower", "level_upper"
        ]
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(df)

    def test_init_inverted_level_numbers_raises(self):
        """
        AssertionError raised when level_number_lower >= level_number_upper.

        The convention is that lower levels have smaller quantum numbers,
        so level_number_lower must always be strictly less than level_number_upper.
        """
        index = pd.MultiIndex.from_tuples(
            [(1, 0, 1, 0)],  # lower=1, upper=0 — inverted
            names=EXPECTED_INDEX_NAMES,
        )
        df = pd.DataFrame(
            {"A_ul": [1e8], "B_ul": [1e-5], "B_lu": [1e-5], "nu": [2.47e15]},
            index=index,
        )
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(df)

    def test_init_equal_level_numbers_raises(self):
        """AssertionError raised when level_number_lower == level_number_upper."""
        index = pd.MultiIndex.from_tuples(
            [(1, 0, 0, 0)],  # lower == upper
            names=EXPECTED_INDEX_NAMES,
        )
        df = pd.DataFrame(
            {"A_ul": [1e8], "B_ul": [1e-5], "B_lu": [1e-5], "nu": [2.47e15]},
            index=index,
        )
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(df)

    def test_init_sorts_index(self, minimal_einstein_coefficients):
        """
        Einstein coefficients are sorted by index after construction.

        Sorted access is required for correct DataFrame alignment when
        computing rates in solve().
        """
        # add a second row with a lower index to force re-sorting
        extra_index = pd.MultiIndex.from_tuples(
            [(1, 0, 2, 3)], names=EXPECTED_INDEX_NAMES
        )
        extra_row = pd.DataFrame(
            {"A_ul": [1e7], "B_ul": [1e-6], "B_lu": [1e-6], "nu": [1e15]},
            index=extra_index,
        )
        unsorted = pd.concat([extra_row, minimal_einstein_coefficients])
        solver = RadiativeRatesSolver(unsorted)
        assert solver.einstein_coefficients.index.is_monotonic_increasing


# ---------------------------------------------------------------------------
# solve() — Output Structure Tests
# ---------------------------------------------------------------------------


class TestRadiativeRatesSolverOutputStructure:
    """Tests for the structure of the DataFrame returned by solve()."""

    def test_solve_returns_dataframe(self, minimal_solver, minimal_rad_field):
        """solve() returns a pandas DataFrame."""
        result = minimal_solver.solve(minimal_rad_field)
        assert isinstance(result, pd.DataFrame)

    def test_solve_output_index_names(self, minimal_solver, minimal_rad_field):
        """
        Output DataFrame has the correct 6-level MultiIndex.

        The index must encode both source and destination ion/level numbers
        so it can be consumed directly by the rate matrix solver.
        """
        result = minimal_solver.solve(minimal_rad_field)
        assert list(result.index.names) == EXPECTED_OUTPUT_INDEX_NAMES

    def test_solve_output_row_count(self, minimal_solver, minimal_rad_field):
        """
        For N transitions, solve() produces 2N rows: one r_lu and one r_ul.

        Both upward (absorption) and downward (emission) transitions must
        appear for each spectral line.
        """
        result = minimal_solver.solve(minimal_rad_field)
        n_transitions = len(minimal_solver.einstein_coefficients)
        assert len(result) == 2 * n_transitions

    def test_solve_output_columns_are_shells(self, minimal_solver, minimal_rad_field):
        """Output columns correspond to the number of shells in the radiation field."""
        result = minimal_solver.solve(minimal_rad_field)
        n_shells = len(minimal_rad_field.temperature)
        assert result.shape[1] == n_shells


# ---------------------------------------------------------------------------
# solve() — Physical Formula Tests
# ---------------------------------------------------------------------------


class TestRadiativeRatesSolverFormulas:
    """Tests verifying the radiative rate equations are correctly implemented."""

    def test_r_lu_equals_B_lu_times_J(
        self, minimal_solver, minimal_rad_field, minimal_einstein_coefficients
    ):
        """
        r_lu = B_lu * J_nu.

        This is the absorption rate (lower to upper). It depends only on
        stimulated absorption — no spontaneous term exists going upward.
        """
        result = minimal_solver.solve(minimal_rad_field)

        B_lu = minimal_einstein_coefficients["B_lu"].values[0]
        nu = minimal_einstein_coefficients["nu"].values[0]
        J_nu = minimal_rad_field.calculate_mean_intensity(
            np.array([nu]) * u.Hz
        )

        expected_r_lu = B_lu * J_nu[0]  # single freq, single shell

        # r_lu row: source=lower(0), dest=upper(1)
        r_lu_row = result.loc[(1, 0, 0, 0, 0, 1)]
        npt.assert_allclose(r_lu_row.values, expected_r_lu, rtol=1e-10)

    def test_r_ul_equals_B_ul_times_J_plus_A_ul(
        self, minimal_solver, minimal_rad_field, minimal_einstein_coefficients
    ):
        """
        r_ul = B_ul * J_nu + A_ul.

        This is the de-excitation rate (upper to lower). It includes both
        stimulated emission (B_ul * J) and spontaneous emission (A_ul).
        The spontaneous term A_ul ensures downward transitions can occur
        even with zero radiation field.
        """
        result = minimal_solver.solve(minimal_rad_field)

        A_ul = minimal_einstein_coefficients["A_ul"].values[0]
        B_ul = minimal_einstein_coefficients["B_ul"].values[0]
        nu = minimal_einstein_coefficients["nu"].values[0]
        J_nu = minimal_rad_field.calculate_mean_intensity(
            np.array([nu]) * u.Hz
        )

        expected_r_ul = B_ul * J_nu[0] + A_ul  # single freq, single shell

        # r_ul row: source=upper(1), dest=lower(0)
        r_ul_row = result.loc[(1, 0, 0, 0, 1, 0)]
        npt.assert_allclose(r_ul_row.values, expected_r_ul, rtol=1e-10)

    def test_zero_radiation_field_r_lu_is_zero(
        self, minimal_einstein_coefficients
    ):
        """
        With W=0, r_lu must be exactly zero.

        There is no radiation to absorb, so no upward transitions can
        occur via radiative absorption.
        """
        rad_field = DilutePlanckianRadiationField(
            temperature=np.array([10000.0]) * u.K,
            dilution_factor=np.array([0.0]),
        )
        solver = RadiativeRatesSolver(minimal_einstein_coefficients)
        result = solver.solve(rad_field)

        r_lu_row = result.loc[(1, 0, 0, 0, 0, 1)]
        npt.assert_array_equal(r_lu_row.values, 0.0)

    def test_zero_radiation_field_r_ul_equals_A_ul(
        self, minimal_einstein_coefficients
    ):
        """
        With W=0, r_ul must equal exactly A_ul (spontaneous emission only).

        Without any radiation field, stimulated emission is zero and only
        spontaneous emission drives downward transitions.
        """
        rad_field = DilutePlanckianRadiationField(
            temperature=np.array([10000.0]) * u.K,
            dilution_factor=np.array([0.0]),
        )
        A_ul = minimal_einstein_coefficients["A_ul"].values[0]
        solver = RadiativeRatesSolver(minimal_einstein_coefficients)
        result = solver.solve(rad_field)

        r_ul_row = result.loc[(1, 0, 0, 0, 1, 0)]
        npt.assert_allclose(r_ul_row.values, A_ul, rtol=1e-10)

    def test_r_ul_always_geq_A_ul(
        self, minimal_solver, minimal_rad_field, minimal_einstein_coefficients
    ):
        """
        r_ul >= A_ul for any non-negative radiation field.

        Since r_ul = B_ul * J + A_ul and both B_ul and J are non-negative,
        the stimulated term can only add to the spontaneous rate.
        """
        result = minimal_solver.solve(minimal_rad_field)
        A_ul = minimal_einstein_coefficients["A_ul"].values[0]
        r_ul_row = result.loc[(1, 0, 0, 0, 1, 0)]
        assert np.all(r_ul_row.values >= A_ul)

    def test_all_rates_non_negative(self, minimal_solver, minimal_rad_field):
        """All rates must be non-negative for a physical radiation field."""
        result = minimal_solver.solve(minimal_rad_field)
        assert np.all(result.values >= 0)

    def test_planckian_radiation_field_accepted(
        self, minimal_solver, minimal_einstein_coefficients
    ):
        """
        solve() works with PlanckianRadiationField (W=1) as input.

        The solver must accept both dilute and pure Planckian radiation fields,
        as they share the same calculate_mean_intensity interface.
        """
        planck_field = PlanckianRadiationField(
            temperature=np.array([10000.0]) * u.K
        )
        result = minimal_solver.solve(planck_field)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2  # one r_lu + one r_ul

    def test_planckian_equals_unit_dilution(
        self, minimal_einstein_coefficients
    ):
        """
        PlanckianRadiationField gives same rates as DilutePlanckian with W=1.

        This validates that the two radiation field types are consistent
        at the W=1 physical limit.
        """
        planck = PlanckianRadiationField(
            temperature=np.array([10000.0]) * u.K
        )
        dilute_unit = DilutePlanckianRadiationField(
            temperature=np.array([10000.0]) * u.K,
            dilution_factor=np.array([1.0]),
        )
        solver = RadiativeRatesSolver(minimal_einstein_coefficients)

        result_planck = solver.solve(planck)
        result_dilute = solver.solve(dilute_unit)

        npt.assert_allclose(
            result_planck.values, result_dilute.values, rtol=1e-10
        )

    def test_rates_scale_with_dilution_factor(
        self, minimal_einstein_coefficients
    ):
        """
        Doubling dilution factor doubles r_lu and increases r_ul by B_ul * delta_J.

        Since r_lu = B_lu * W * B_nu(T), it scales linearly with W.
        r_ul = B_ul * W * B_nu(T) + A_ul, so the stimulated part scales linearly
        but the total does not double due to the constant A_ul offset.
        """
        rad_field_1 = DilutePlanckianRadiationField(
            temperature=np.array([10000.0]) * u.K,
            dilution_factor=np.array([0.2]),
        )
        rad_field_2 = DilutePlanckianRadiationField(
            temperature=np.array([10000.0]) * u.K,
            dilution_factor=np.array([0.4]),  # doubled
        )
        solver = RadiativeRatesSolver(minimal_einstein_coefficients)

        result_1 = solver.solve(rad_field_1)
        result_2 = solver.solve(rad_field_2)

        # r_lu doubles exactly
        r_lu_1 = result_1.loc[(1, 0, 0, 0, 0, 1)].values
        r_lu_2 = result_2.loc[(1, 0, 0, 0, 0, 1)].values
        npt.assert_allclose(r_lu_2, 2 * r_lu_1, rtol=1e-10)


# ---------------------------------------------------------------------------
# solve() — Regression Tests (using real Si Chianti atomic data)
# ---------------------------------------------------------------------------


class TestRadiativeRatesSolverRegression:
    """Regression tests pinning solve() output against stored reference data."""

    def test_solve_regression_dilute_planckian(
        self, radiative_rate_solver, regression_data
    ):
        """
        Pin full solve() output for Si I lines against regression data.

        Uses a DilutePlanckianRadiationField with realistic supernova
        photosphere conditions (T=10000 K, W=0.5) over 20 shells.
        Any change to the rate equations, index construction, or unit
        handling will be caught immediately by this test.
        """
        n_shells = 20
        rad_field = DilutePlanckianRadiationField(
            temperature=np.ones(n_shells) * 10000.0 * u.K,
            dilution_factor=np.ones(n_shells) * 0.5,
        )
        actual = radiative_rate_solver.solve(rad_field)
        expected = regression_data.sync_dataframe(
            actual, key="radiative_rates_dilute_planckian"
        )
        pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-13)

    def test_solve_regression_pure_planckian(
        self, radiative_rate_solver, regression_data
    ):
        """
        Pin solve() output for a pure PlanckianRadiationField (W=1).

        This tests the W=1 limit separately from the dilute case
        to confirm both radiation field types produce consistent,
        reproducible results.
        """
        n_shells = 20
        rad_field = PlanckianRadiationField(
            temperature=np.ones(n_shells) * 10000.0 * u.K
        )
        actual = radiative_rate_solver.solve(rad_field)
        expected = regression_data.sync_dataframe(
            actual, key="radiative_rates_pure_planckian"
        )
        pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-13)
