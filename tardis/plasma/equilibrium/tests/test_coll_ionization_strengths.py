import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)


@pytest.fixture
def collisional_ionization_solver(mock_photoionization_cross_sections):
    """Fixture for CollisionalIonizationSeaton instance."""
    return CollisionalIonizationSeaton(mock_photoionization_cross_sections)


def test_solve(collisional_ionization_solver):
    """Test the solve method with a valid electron temperature."""
    electron_temperature = 1e4 * u.K
    result = collisional_ionization_solver.solve(electron_temperature)

    assert isinstance(result, pd.DataFrame)
    assert not result.empty
    assert result.index.equals(
        collisional_ionization_solver.photoionization_cross_sections.index
    )
    assert np.all(result.values >= 0)
