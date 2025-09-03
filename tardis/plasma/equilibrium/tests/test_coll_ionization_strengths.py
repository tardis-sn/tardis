import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)


@pytest.fixture
def collisional_ionization_solver(mock_photoionization_cross_sections):
    """Fixture for CollisionalIonizationSeaton instance."""
    return CollisionalIonizationSeaton(mock_photoionization_cross_sections)


def test_solve(collisional_ionization_solver, regression_data):
    """Test the solve method with a valid electron temperature."""
    electron_temperature = 1e4 * u.K
    actual_result = collisional_ionization_solver.solve(electron_temperature)

    assert isinstance(actual_result, pd.DataFrame)
    assert not actual_result.empty
    assert actual_result.index.equals(
        collisional_ionization_solver.photoionization_cross_sections.index
    )
    assert np.all(actual_result.values >= 0)

    # Regression data comparison
    expected_result = regression_data.sync_dataframe(
        actual_result, key="collisional_ionization_seaton_solve"
    )

    pdt.assert_frame_equal(actual_result, expected_result)
