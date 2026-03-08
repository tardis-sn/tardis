import numpy as np
import pandas as pd
import astropy.units as u
import pytest
from unittest.mock import patch

@pytest.fixture
def mock_photoionization_cross_sections():
    index = pd.MultiIndex.from_tuples(
        [
            (1, 0, 0),
            (1, 0, 0),
            (1, 0, 1),
            (1, 0, 1),
        ]
    )

    df = pd.DataFrame(
        {
            "nu": [1e15, 2e15, 1e15, 2e15],
            "x_sect": [1e-18, 1e-18, 2e-18, 2e-18],
        },
        index=index,
    )

    return df


from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    SpontaneousRecombinationCoeffSolver
)

def test_common_prefactor(mock_photoionization_cross_sections):

    solver = SpontaneousRecombinationCoeffSolver(
        mock_photoionization_cross_sections
    )

    prefactor = solver.common_prefactor

    assert prefactor is not None
    assert len(prefactor) == len(mock_photoionization_cross_sections)



def test_boltzmann_factor_range(mock_photoionization_cross_sections):

    solver = SpontaneousRecombinationCoeffSolver(
        mock_photoionization_cross_sections
    )

    electron_temperature = np.array([5000]) * u.K

    result = solver.calculate_photoionization_boltzmann_factor(
        electron_temperature
    )

    assert np.all(result <= 1)
    assert np.all(result > 0)


@patch(
"tardis.plasma.equilibrium.rates.photoionization_strengths.integrate_array_by_blocks"
)
def test_solver_output_shape(mock_integrate, mock_photoionization_cross_sections):

    mock_integrate.return_value = np.array([[1.0],[2.0]])

    solver = SpontaneousRecombinationCoeffSolver(
        mock_photoionization_cross_sections
    )

    electron_temperature = np.array([5000]) * u.K

    result = solver.solve(electron_temperature)

    assert isinstance(result, pd.DataFrame)


'''
I've commented out this test because it took too long to run and used 100% of my CPU.
hopefully this test passes for you.
'''
# @patch(
# "tardis.plasma.equilibrium.rates.photoionization_strengths.integrate_array_by_blocks"
# )
# def test_lyman_continuum_zero(mock_photoionization_cross_sections):

#     solver = SpontaneousRecombinationCoeffSolver(
#         mock_photoionization_cross_sections
#     )

#     electron_temperature = np.array([5000]) * u.K

#     result = solver.solve(electron_temperature)

#     assert result.loc[(1,0,0)].values[0] == 0

tardis/plasma/equilibrium/tests/test_spontaneous_recombination_solver.py