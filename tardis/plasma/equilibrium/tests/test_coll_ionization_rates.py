import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.equilibrium.rates.collisional_ionization_rates import (
    CollisionalIonizationRateSolver,
)


class MockElectronDistribution:
    """Mock class for electron distribution."""

    def __init__(self, temperature, number_density):
        self.temperature = temperature
        self.number_density = number_density


@pytest.fixture
def mock_electron_distribution():
    """Fixture for mock electron distribution."""
    return MockElectronDistribution(
        temperature=np.array([10000, 20000]) * u.K,
        number_density=np.array([1e13, 2e13]) * u.cm**-3,
    )


@pytest.fixture
def mock_saha_factor():
    """Fixture for mock Saha factor."""
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    return pd.DataFrame([[1e-15, 1e-20], [1e-15, 1e-20]], index=index)


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
    regression_data,
):
    """Test the solve method of CollisionalIonizationRateSolver."""
    solver = CollisionalIonizationRateSolver(
        mock_photoionization_cross_sections
    )

    actual_ionization_rates, actual_recombination_rates = solver.solve(
        mock_electron_distribution, mock_saha_factor, approximation="seaton"
    )

    # write paths manually with regression_data directory info from the class
    if regression_data.enable_generate_reference:
        actual_ionization_rates.to_hdf(
            regression_data.absolute_regression_data_dir
            / "ionization_rates.h5",
            key="data",
        )
        actual_recombination_rates.to_hdf(
            regression_data.absolute_regression_data_dir
            / "recombination_rates.h5",
            key="data",
        )
        pytest.skip("Skipping test to generate reference data")
    else:
        expected_ionization_rates = pd.read_hdf(
            regression_data.absolute_regression_data_dir
            / "ionization_rates.h5",
            key="data",
        )

        expected_recombination_rates = pd.read_hdf(
            regression_data.absolute_regression_data_dir
            / "recombination_rates.h5",
            key="data",
        )

        pdt.assert_frame_equal(
            actual_ionization_rates, expected_ionization_rates
        )
        pdt.assert_frame_equal(
            actual_recombination_rates, expected_recombination_rates
        )


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
