import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)


@pytest.fixture
def mock_photoionization_cross_sections():
    """Fixture for mock photoionization cross-sections."""
    data = {
        "nu": [1e15, 2e15, 3e15],
        "x_sect": [1e-18, 2e-18, 3e-18],
    }
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 1, 0), (2, 0, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    return pd.DataFrame(data, index=index)


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
