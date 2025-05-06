import pytest
import pandas as pd
import numpy as np

from tardis.plasma.equilibrium.rates.collisional_ionization_rates import (
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.collisional_ionization_strengths import (
    CollisionalIonizationSeaton,
)


class MockElectronDistribution:
    """Mock class for electron distribution."""

    def __init__(self, temperature, number_density):
        self.temperature = temperature
        self.number_density = number_density


@pytest.fixture
def mock_photoionization_cross_sections():
    """Fixture for mock photoionization cross-sections."""
    return pd.DataFrame(
        {
            "atomic_number": [1, 2],
            "ion_number": [0, 1],
            "cross_section": [1e-18, 2e-18],
        }
    )


@pytest.fixture
def mock_electron_distribution():
    """Fixture for mock electron distribution."""
    return MockElectronDistribution(
        temperature=pd.Series([10000, 20000], index=[0, 1]),
        number_density=pd.Series([1e13, 2e13], index=[0, 1]),
    )


@pytest.fixture
def mock_saha_factor():
    """Fixture for mock Saha factor."""
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (2, 1, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    return pd.DataFrame([0.1, 0.2], index=index, columns=["saha_factor"])


def test_collisional_ionization_rate_solver_init(
    mock_photoionization_cross_sections,
):
    """Test initialization of CollisionalIonizationRateSolver."""
    solver = CollisionalIonizationRateSolver(
        mock_photoionization_cross_sections
    )
    assert solver.photoionization_cross_sections.equals(
        mock_photoionization_cross_sections
    )


def test_collisional_ionization_rate_solver_solve(
    mock_photoionization_cross_sections,
    mock_electron_distribution,
    mock_saha_factor,
):
    """Test the solve method of CollisionalIonizationRateSolver."""
    solver = CollisionalIonizationRateSolver(
        mock_photoionization_cross_sections
    )

    # Mock the strength solver
    class MockStrengthSolver:
        def __init__(self, photoionization_cross_sections):
            pass

        def solve(self, temperature):
            return pd.DataFrame(
                {
                    "atomic_number": [1, 2],
                    "ion_number": [0, 1],
                    "level_number_source": [0, 0],
                    "rate": [1e-8, 2e-8],
                }
            ).set_index(["atomic_number", "ion_number", "level_number_source"])

    # Replace the strength solver with the mock
    CollisionalIonizationSeaton.solve = MockStrengthSolver.solve

    ionization_rates, recombination_rates = solver.solve(
        mock_electron_distribution, mock_saha_factor, approximation="seaton"
    )

    # Check ionization rates
    assert not ionization_rates.empty
    assert "rate" in ionization_rates.columns

    # Check recombination rates
    assert not recombination_rates.empty
    assert "rate" in recombination_rates.columns


def test_collisional_ionization_rate_solver_invalid_approximation(
    mock_photoionization_cross_sections,
    mock_electron_distribution,
    mock_saha_factor,
):
    """Test that an invalid approximation raises a ValueError."""
    solver = CollisionalIonizationRateSolver(
        mock_photoionization_cross_sections
    )
    with pytest.raises(
        ValueError, match="approximation invalid_approx not supported"
    ):
        solver.solve(
            mock_electron_distribution,
            mock_saha_factor,
            approximation="invalid_approx",
        )
