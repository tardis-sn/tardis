import astropy.units as u
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.io.atom_data import AtomData
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

@pytest.mark.xfail(strict=True, raises=AssertionError)
def test_invalid_coefficients(invalid_coefficients):
    solver = RadiativeRatesSolver(invalid_coefficients)
    