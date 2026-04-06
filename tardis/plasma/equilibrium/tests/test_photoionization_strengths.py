import astropy.units as u
import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    AnalyticCorrectedPhotoionizationCoeffSolver,
    AnalyticPhotoionizationCoeffSolver,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
)
from tardis.plasma.radiation_field import (
    DilutePlanckianRadiationField,
)
from tardis.transport.montecarlo.estimators import init_estimators_continuum


@pytest.fixture
def mock_dilute_blackbody_radiationfield_state():
    return DilutePlanckianRadiationField(
        np.ones(20) * 10000 * u.K, dilution_factor=np.ones(20) * 0.5
    )


def test_spontaneous_recombination_coeff_solver(
    mock_photoionization_cross_sections,
    regression_data,
):
    solver = SpontaneousRecombinationCoeffSolver(
        mock_photoionization_cross_sections
    )
    electron_temperature = np.array([1e4, 1e4]) * u.K
    actual_result = solver.solve(electron_temperature)

    assert isinstance(actual_result, pd.DataFrame)
    assert actual_result.shape[0] == len(
        mock_photoionization_cross_sections.index.unique()
    )
    assert not actual_result.isnull().values.any()

    # Regression data comparison
    expected_result = regression_data.sync_dataframe(
        actual_result, key="spontaneous_recombination_coeff"
    )

    pdt.assert_frame_equal(actual_result, expected_result, atol=0, rtol=1e-15)


def test_analytic_photoionization_coeff_solver(
    mock_photoionization_cross_sections,
    mock_dilute_blackbody_radiationfield_state,
    regression_data,
):
    solver = AnalyticPhotoionizationCoeffSolver(
        mock_photoionization_cross_sections
    )
    electron_temperature = np.array([1e4] * 20) * u.K

    (
        actual_photoionization_rate_coeff,
        actual_stimulated_recombination_rate_coeff,
    ) = solver.solve(
        mock_dilute_blackbody_radiationfield_state, electron_temperature
    )

    assert isinstance(actual_photoionization_rate_coeff, pd.DataFrame)
    assert isinstance(actual_stimulated_recombination_rate_coeff, pd.DataFrame)
    assert actual_photoionization_rate_coeff.shape[0] == len(
        mock_photoionization_cross_sections.index.unique()
    )
    assert actual_stimulated_recombination_rate_coeff.shape[0] == len(
        mock_photoionization_cross_sections.index.unique()
    )

    # Regression data comparison
    expected_photoionization_rate_coeff = regression_data.sync_dataframe(
        actual_photoionization_rate_coeff,
        key="analytic_photoionization_rate_coeff",
    )
    expected_stimulated_recombination_rate_coeff = (
        regression_data.sync_dataframe(
            actual_stimulated_recombination_rate_coeff,
            key="analytic_stimulated_recombination_rate_coeff",
        )
    )

    pdt.assert_frame_equal(
        actual_photoionization_rate_coeff,
        expected_photoionization_rate_coeff,
        atol=0,
        rtol=1e-15,
    )
    pdt.assert_frame_equal(
        actual_stimulated_recombination_rate_coeff,
        expected_stimulated_recombination_rate_coeff,
        atol=0,
        rtol=1e-15,
    )


def test_analytic_corrected_photoionization_coeff_solver(
    mock_photoionization_cross_sections,
    mock_dilute_blackbody_radiationfield_state,
    regression_data,
):
    solver = AnalyticCorrectedPhotoionizationCoeffSolver(
        mock_photoionization_cross_sections
    )
    electron_temperature = np.array([1e4, 1e4]) * u.K
    lte_level_population = pd.DataFrame(
        np.ones((2, 2)), index=mock_photoionization_cross_sections.index
    )
    level_population = pd.DataFrame(
        np.ones((2, 2)), index=mock_photoionization_cross_sections.index
    )
    lte_ion_population = pd.DataFrame(
        np.ones((2, 2)), index=mock_photoionization_cross_sections.index
    )
    ion_population = pd.DataFrame(
        np.ones((2, 2)), index=mock_photoionization_cross_sections.index
    )

    actual_corrected_photoionization_rate_coeff = solver.solve(
        mock_dilute_blackbody_radiationfield_state,
        electron_temperature,
        lte_level_population,
        level_population,
        lte_ion_population,
        ion_population,
    )

    assert isinstance(actual_corrected_photoionization_rate_coeff, pd.DataFrame)
    assert actual_corrected_photoionization_rate_coeff.shape[0] == len(
        mock_photoionization_cross_sections.index.unique()
    )

    # Regression data comparison
    expected_corrected_photoionization_rate_coeff = (
        regression_data.sync_dataframe(
            actual_corrected_photoionization_rate_coeff,
            key="analytic_corrected_photoionization_rate_coeff",
        )
    )

    pdt.assert_frame_equal(
        actual_corrected_photoionization_rate_coeff,
        expected_corrected_photoionization_rate_coeff,
        atol=0,
        rtol=1e-15,
    )


def test_estimated_photoionization_coeff_solver(regression_data):
    level2continuum_edge_idx = pd.Series(np.array([0, 2, 4]))
    solver = EstimatedPhotoionizationCoeffSolver(level2continuum_edge_idx)

    # Create mock estimator arrays
    photo_ion_array = np.array([[1e-5], [2e-5], [3e-5]])
    stim_recomb_array = np.array([[1e-6], [2e-6], [3e-6]])

    time_simulation = 1e5 * u.s
    volume = 1e30 * u.cm**3

    # Initialize with correct factory signature
    estimators_continuum = init_estimators_continuum(
        n_levels_bf_species_by_n_cells_tuple=(3, 1), n_cells=1
    )

    # Set the estimator values
    estimators_continuum.photo_ion_estimator[:] = photo_ion_array
    estimators_continuum.stim_recomb_estimator[:] = stim_recomb_array

    (
        actual_photoionization_rate_coeff,
        actual_stimulated_recombination_rate_coeff,
    ) = solver.solve(estimators_continuum, time_simulation, volume)

    assert isinstance(actual_photoionization_rate_coeff, pd.DataFrame)
    assert isinstance(actual_stimulated_recombination_rate_coeff, pd.DataFrame)
    assert actual_photoionization_rate_coeff.shape[0] == len(
        level2continuum_edge_idx
    )
    assert actual_stimulated_recombination_rate_coeff.shape[0] == len(
        level2continuum_edge_idx
    )

    # Regression data comparison
    expected_photoionization_rate_coeff = regression_data.sync_dataframe(
        actual_photoionization_rate_coeff,
        key="estimated_photoionization_rate_coeff",
    )
    expected_stimulated_recombination_rate_coeff = (
        regression_data.sync_dataframe(
            actual_stimulated_recombination_rate_coeff,
            key="estimated_stimulated_recombination_rate_coeff",
        )
    )

    pdt.assert_frame_equal(
        actual_photoionization_rate_coeff,
        expected_photoionization_rate_coeff,
        atol=0,
        rtol=1e-15,
    )
    pdt.assert_frame_equal(
        actual_stimulated_recombination_rate_coeff,
        expected_stimulated_recombination_rate_coeff,
        atol=0,
        rtol=1e-15,
    )
