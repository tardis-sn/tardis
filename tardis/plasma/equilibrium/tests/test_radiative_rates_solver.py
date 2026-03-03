import numpy as np
import pandas as pd
import pytest

from tardis.plasma.equilibrium.rates.radiative_rates import (
    RadiativeRatesSolver,
)


def make_valid_einstein_df():
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    )

    return pd.DataFrame(
        {
            "A_ul": [1.0],
            "B_ul": [2.0],
            "B_lu": [3.0],
            "nu": [4.0],
        },
        index=index,
    )


def test_missing_required_column_raises():
    df = make_valid_einstein_df().drop(columns=["A_ul"])

    with pytest.raises(AssertionError):
        RadiativeRatesSolver(df)


class DummyRadiationField:
    def calculate_mean_intensity(self, nu):
        return np.ones_like(nu)


def test_solve_returns_dataframe():
    df = make_valid_einstein_df()
    solver = RadiativeRatesSolver(df)

    rad_field = DummyRadiationField()
    result = solver.solve(rad_field)

    assert isinstance(result, pd.DataFrame)
    assert not result.empty


def test_solve_output_index_structure():
    df = make_valid_einstein_df()
    solver = RadiativeRatesSolver(df)

    rad_field = DummyRadiationField()
    result = solver.solve(rad_field)

    expected_index_names = [
        "atomic_number",
        "ion_number",
        "ion_number_source",
        "ion_number_destination",
        "level_number_source",
        "level_number_destination",
    ]

    assert list(result.index.names) == expected_index_names