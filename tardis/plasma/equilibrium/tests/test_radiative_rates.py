import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.equilibrium.rates.radiative_rates import RadiativeRatesSolver
from tardis.plasma.radiation_field.planck_rad_field import (
    PlanckianRadiationField,
)

invalid_index_df = pd.DataFrame(
    {
        "A_ul": [1e8],
        "B_ul": [1e-19],
        "B_lu": [2e-19],
        "nu": [3e15],
    },
    index=pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=[
            "atomic_num",
            "ion_numb",
            "level_lower",
            "level_upper",
        ],
    ),
)

invalid_column_df = pd.DataFrame(
    {
        "A_ul": [1e8],
        "B_ul": [1e-19],
        "B_lu": [2e-19],
    },
    index=pd.MultiIndex.from_tuples(
        [(1, 0, 0, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    ),
)

invalid_lower_higher_df = pd.DataFrame(
    {
        "A_ul": [1e8],
        "B_ul": [1e-19],
        "B_lu": [2e-19],
        "nu": [3e15],
    },
    index=pd.MultiIndex.from_tuples(
        [(1, 0, 2, 1)],
        names=[
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ],
    ),
)


@pytest.fixture(
    scope="function",
    params=[invalid_index_df, invalid_column_df, invalid_lower_higher_df]
)
def invalid_coefficients(request):
    return request.param


@pytest.fixture(scope="class")
def mock_radiation_field():
    """Fixture for mock radiation field."""
    temperature = [10000, 20000]* u.K
    return PlanckianRadiationField(temperature=temperature)

def test_radiative_rate_solver_init(new_chianti_atomic_dataset,regression_data):
    einstein_coefficients_df = new_chianti_atomic_dataset.lines.xs((1,0),drop_level=False)
    solver = RadiativeRatesSolver(einstein_coefficients_df)
    actual_einstein_coeffs = solver.einstein_coefficients
    expected_einstein_coeffs = regression_data.sync_dataframe(
        actual_einstein_coeffs, key="einstein_coeffs")
    # CAN DO A NORMAL ASSERT WHEN ATOMIC DATA IS UPDATED FOR PANDAS 3.X
    assert actual_einstein_coeffs.columns.names == ["N."]
    pdt.assert_frame_equal(
        actual_einstein_coeffs, expected_einstein_coeffs, check_names=False
    )

def test_radiative_rate_solver_solve(new_chianti_atomic_dataset, mock_radiation_field, regression_data):
    einstein_coefficients_df = new_chianti_atomic_dataset.lines.xs((1,0),drop_level=False)
    solver = RadiativeRatesSolver(einstein_coefficients_df)
    actual_radiative_rates = solver.solve(mock_radiation_field)
    expected_radiative_rates = regression_data.sync_dataframe(
        actual_radiative_rates, key="solved_radiative_rates"
    )
    pdt.assert_frame_equal(actual_radiative_rates,expected_radiative_rates,atol=0,rtol=1e-15)


def test_radiative_rate_solver_applies_sobolev_escape_probability(
    new_chianti_atomic_dataset,
):
    einstein_coefficients = new_chianti_atomic_dataset.lines.loc[
        (1, 0, slice(None), slice(None)),
        ["A_ul", "B_ul", "B_lu", "nu"],
    ].iloc[:1]
    solver = RadiativeRatesSolver(einstein_coefficients)
    radiation_field = PlanckianRadiationField(
        temperature=np.array([10000.0, 10000.0]) * u.K
    )
    beta_sobolevs = pd.DataFrame(
        [[0.25, 0.5]], index=einstein_coefficients.index, columns=[0, 1]
    )

    unscaled_rates = solver.solve(radiation_field)
    scaled_rates = solver.solve(radiation_field, beta_sobolevs=beta_sobolevs)

    pdt.assert_frame_equal(
        scaled_rates,
        unscaled_rates.multiply(beta_sobolevs.iloc[0].to_numpy(), axis=1),
    )

@pytest.mark.xfail(strict=True, raises=AssertionError)
def test_invalid_coefficients(invalid_coefficients):
    solver = RadiativeRatesSolver(invalid_coefficients)

