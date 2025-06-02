import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticCorrectedPhotoionizationCoeffSolver,
    AnalyticPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)


@pytest.fixture
def mock_photoionization_cross_sections():
    """Fixture to create mock photoionization cross-sections."""
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1), (2, 1, 0), (2, 1, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    data = {
        "nu": [1e15, 2e15, 1e15, 2e15],
        "x_sect": [1e-18, 2e-18, 1e-18, 2e-18],
    }
    return pd.DataFrame(data, index=index)


@pytest.fixture
def mock_dilute_blackbody_radiationfield_state():
    """Fixture to create a mock radiation field state."""

    class MockRadiationFieldState:
        def __init__(self):
            self.temperature = np.array([1e4, 1e4]) * u.K

        def calculate_mean_intensity(self, nu):
            return np.array([1e-5, 2e-5]) * u.erg / (u.cm**2 * u.s * u.Hz)

    return MockRadiationFieldState()


def test_spontaneous_recombination_coeff_solver(
    mock_photoionization_cross_sections,
):
    solver = SpontaneousRecombinationCoeffSolver(
        mock_photoionization_cross_sections
    )
    electron_temperature = np.array([1e4, 1e4]) * u.K
    result = solver.solve(electron_temperature)

    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == len(
        mock_photoionization_cross_sections.index.unique()
    )
    assert not result.isnull().values.any()


def test_analytic_photoionization_coeff_solver(
    mock_photoionization_cross_sections,
    mock_dilute_blackbody_radiationfield_state,
):
    solver = AnalyticPhotoionizationCoeffSolver(
        mock_photoionization_cross_sections
    )
    electron_temperature = np.array([1e4, 1e4]) * u.K

    photoionization_rate_coeff, stimulated_recombination_rate_coeff = (
        solver.solve(
            mock_dilute_blackbody_radiationfield_state, electron_temperature
        )
    )

    assert isinstance(photoionization_rate_coeff, pd.DataFrame)
    assert isinstance(stimulated_recombination_rate_coeff, pd.DataFrame)
    assert photoionization_rate_coeff.shape[0] == len(
        mock_photoionization_cross_sections.index.unique()
    )
    assert stimulated_recombination_rate_coeff.shape[0] == len(
        mock_photoionization_cross_sections.index.unique()
    )


def test_analytic_corrected_photoionization_coeff_solver(
    mock_photoionization_cross_sections,
    mock_dilute_blackbody_radiationfield_state,
):
    solver = AnalyticCorrectedPhotoionizationCoeffSolver(
        mock_photoionization_cross_sections
    )
    electron_temperature = np.array([1e4, 1e4]) * u.K
    lte_level_population = pd.DataFrame(
        np.ones((4, 2)), index=mock_photoionization_cross_sections.index
    )
    level_population = pd.DataFrame(
        np.ones((4, 2)), index=mock_photoionization_cross_sections.index
    )
    lte_ion_population = pd.DataFrame(
        np.ones((4, 2)), index=mock_photoionization_cross_sections.index
    )
    ion_population = pd.DataFrame(
        np.ones((4, 2)), index=mock_photoionization_cross_sections.index
    )

    corrected_photoionization_rate_coeff = solver.solve(
        mock_dilute_blackbody_radiationfield_state,
        electron_temperature,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
    )

    assert isinstance(corrected_photoionization_rate_coeff, pd.DataFrame)
    assert corrected_photoionization_rate_coeff.shape[0] == len(
        mock_photoionization_cross_sections.index.unique()
    )


def test_estimated_photoionization_coeff_solver():
    level2continuum_edge_idx = np.array([0, 2, 4])
    solver = EstimatedPhotoionizationCoeffSolver(level2continuum_edge_idx)

    class MockRadFieldMCEstimators:
        photo_ion_estimator = np.array([1e-5, 2e-5, 3e-5])
        stim_recomb_estimator = np.array([1e-6, 2e-6, 3e-6])

    radfield_mc_estimators = MockRadFieldMCEstimators()
    time_simulation = 1e5 * u.s
    volume = 1e30 * u.cm**3

    photoionization_rate_coeff, stimulated_recombination_rate_coeff = (
        solver.solve(radfield_mc_estimators, time_simulation, volume)
    )

    assert isinstance(photoionization_rate_coeff, pd.DataFrame)
    assert isinstance(stimulated_recombination_rate_coeff, pd.DataFrame)
    assert (
        photoionization_rate_coeff.shape[0] == len(level2continuum_edge_idx) - 1
    )
    assert (
        stimulated_recombination_rate_coeff.shape[0]
        == len(level2continuum_edge_idx) - 1
    )
