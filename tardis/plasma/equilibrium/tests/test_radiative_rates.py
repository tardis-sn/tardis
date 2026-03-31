"""
Unit tests for tardis.plasma.equilibrium.rates.radiative_rates.RadiativeRatesSolver.

These tests cover:
- Initialization validation (correct MultiIndex names, required columns,
  level_number ordering)
- solve() output shape and index names
- Physics correctness: with zero dilution factor J_nu = 0, so r_lu = 0 and
  r_ul = A_ul (spontaneous emission only)
- Regression comparison at rtol=1e-15 against stored golden data
- Non-negativity of all computed rates
"""

import numpy as np
import astropy.units as u
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.plasma.equilibrium.rates import RadiativeRatesSolver


# ---------------------------------------------------------------------------
# Helper fixtures local to this module
# ---------------------------------------------------------------------------


@pytest.fixture
def rad_field(collisional_simulation_state):
    """DilutePlanckianRadiationField with W=0.5 in all shells."""
    return DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.full(
            len(collisional_simulation_state.t_radiative), 0.5
        ),
    )


@pytest.fixture
def zero_dilution_rad_field(collisional_simulation_state):
    """DilutePlanckianRadiationField with W=0 (optically thin / dark limit)."""
    return DilutePlanckianRadiationField(
        collisional_simulation_state.t_radiative,
        dilution_factor=np.zeros(
            len(collisional_simulation_state.t_radiative)
        ),
    )


# ---------------------------------------------------------------------------
# Initialization tests
# ---------------------------------------------------------------------------


class TestRadiativeRatesSolverInit:
    """Tests for RadiativeRatesSolver.__init__ validation logic."""

    def test_valid_init_stores_einstein_coefficients(
        self, radiative_rate_solver
    ):
        """Solver should store the einstein_coefficients DataFrame."""
        assert hasattr(radiative_rate_solver, "einstein_coefficients")
        assert {"A_ul", "B_ul", "B_lu", "nu"}.issubset(
            radiative_rate_solver.einstein_coefficients.columns
        )

    def test_wrong_index_names_raises(self):
        """Solver must reject DataFrames whose MultiIndex has wrong level names."""
        bad_df = pd.DataFrame(
            {"A_ul": [1.0], "B_ul": [1.0], "B_lu": [1.0], "nu": [1e14]},
            index=pd.MultiIndex.from_tuples(
                [(1, 0, 0, 1)], names=["a", "b", "c", "d"]
            ),
        )
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(bad_df)

    def test_missing_column_raises(self, radiative_transitions):
        """Solver must reject DataFrames missing any of A_ul/B_ul/B_lu/nu."""
        for col in ["A_ul", "B_ul", "B_lu", "nu"]:
            incomplete = radiative_transitions.drop(columns=[col])
            with pytest.raises(AssertionError):
                RadiativeRatesSolver(incomplete)

    def test_inverted_level_numbers_raises(self, radiative_transitions):
        """Solver must reject transitions where level_number_lower >= level_number_upper."""
        inverted = radiative_transitions.copy()
        inverted.index = inverted.index.swaplevel(
            "level_number_lower", "level_number_upper"
        )
        with pytest.raises(AssertionError):
            RadiativeRatesSolver(inverted)


# ---------------------------------------------------------------------------
# solve() structural tests
# ---------------------------------------------------------------------------


class TestRadiativeRatesSolverSolveStructure:
    """Tests for output structure of RadiativeRatesSolver.solve()."""

    def test_solve_output_shape(self, radiative_rate_solver, rad_field):
        """
        Output must have exactly 2 * N_transitions rows (r_lu + r_ul for each
        transition) and N_shells columns.
        """
        result = radiative_rate_solver.solve(rad_field)
        n_transitions = len(radiative_rate_solver.einstein_coefficients)
        n_shells = len(rad_field.temperature)
        assert result.shape == (2 * n_transitions, n_shells)

    def test_solve_output_index_names(self, radiative_rate_solver, rad_field):
        """Output MultiIndex must have the six expected level names."""
        result = radiative_rate_solver.solve(rad_field)
        expected_names = [
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ]
        assert list(result.index.names) == expected_names

    def test_solve_returns_dataframe(self, radiative_rate_solver, rad_field):
        """solve() must return a pandas DataFrame."""
        result = radiative_rate_solver.solve(rad_field)
        assert isinstance(result, pd.DataFrame)


# ---------------------------------------------------------------------------
# Physics correctness tests
# ---------------------------------------------------------------------------


class TestRadiativeRatesSolverPhysics:
    """
    Physics-verifying tests for RadiativeRatesSolver.solve().

    The radiative transition rates are:
        r_lu = B_lu * J_nu             (absorption)
        r_ul = B_ul * J_nu + A_ul      (stimulated + spontaneous emission)

    where J_nu = W * B_nu(T, nu) is the dilute mean intensity.
    """

    def test_zero_dilution_gives_only_spontaneous_emission(
        self, radiative_rate_solver, zero_dilution_rad_field
    ):
        """
        With W=0: J_nu = 0, so:
          r_lu = 0  (no absorption without radiation)
          r_ul = A_ul  (only spontaneous emission remains)

        This directly validates the physics formula against the Einstein A coefficient.
        """
        result = radiative_rate_solver.solve(zero_dilution_rad_field)
        ec = radiative_rate_solver.einstein_coefficients
        n_shells = len(zero_dilution_rad_field.temperature)

        # Upward transitions (source level < destination level): r_lu must be 0
        r_lu = result.loc[
            result.index.get_level_values("level_number_source")
            < result.index.get_level_values("level_number_destination")
        ]
        np.testing.assert_allclose(
            r_lu.values,
            np.zeros((len(ec), n_shells)),
            atol=1e-30,
            err_msg="Absorption rates must be zero when dilution_factor=0",
        )

        # Downward transitions (source level > destination level): r_ul = A_ul
        r_ul = result.loc[
            result.index.get_level_values("level_number_source")
            > result.index.get_level_values("level_number_destination")
        ]
        expected_r_ul = np.outer(ec["A_ul"].values, np.ones(n_shells))
        np.testing.assert_allclose(
            r_ul.values,
            expected_r_ul,
            rtol=1e-10,
            err_msg="Emission rates must equal A_ul when dilution_factor=0",
        )

    def test_rates_are_non_negative(self, radiative_rate_solver, rad_field):
        """All transition rates must be non-negative (physically required)."""
        result = radiative_rate_solver.solve(rad_field)
        assert (result.values >= 0).all(), (
            "Negative transition rates detected — physical constraint violated."
        )

    def test_higher_dilution_gives_higher_rates(
        self, radiative_rate_solver, collisional_simulation_state
    ):
        """
        Increasing the dilution factor must increase all rates (since J_nu
        scales linearly with W, and A_ul >= 0 ensures r_ul is also increasing).
        """
        t_rad = collisional_simulation_state.t_radiative

        low_W_field = DilutePlanckianRadiationField(
            t_rad, dilution_factor=np.full(len(t_rad), 0.1)
        )
        high_W_field = DilutePlanckianRadiationField(
            t_rad, dilution_factor=np.full(len(t_rad), 0.5)
        )

        low_result = radiative_rate_solver.solve(low_W_field)
        high_result = radiative_rate_solver.solve(high_W_field)

        assert (high_result.values >= low_result.values).all(), (
            "Higher dilution factor must produce higher or equal transition rates."
        )


# ---------------------------------------------------------------------------
# Regression test
# ---------------------------------------------------------------------------


class TestRadiativeRatesSolverRegression:
    """Regression test: compare solve() output against stored golden data."""

    def test_solve_regression(
        self, radiative_rate_solver, rad_field, regression_data
    ):
        """
        Full regression test: solve() output must match stored golden data
        at rtol=1e-15 (machine-precision level).
        """
        actual = radiative_rate_solver.solve(rad_field)
        expected = regression_data.sync_dataframe(actual)
        pdt.assert_frame_equal(actual, expected, atol=0, rtol=1e-15)
