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
    apply_lte_level_to_ion_factor,
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


def test_estimated_photoionization_coeff_solver():
    photoionization_index = pd.MultiIndex.from_tuples(
        [(1, 0, 0), (1, 0, 1), (1, 1, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    level2continuum_edge_idx = pd.Series(
        np.array([0, 2, 4]), index=photoionization_index
    )
    solver = EstimatedPhotoionizationCoeffSolver(level2continuum_edge_idx)

    # Create mock estimator arrays
    photo_ion_array = np.array([[1e-5], [2e-5], [3e-5]])
    stim_recomb_array = np.array([[1e-6], [2e-6], [3e-6]])

    estimator_state = init_estimators_continuum(
        n_levels_bf_species_by_n_cells_tuple=(3, 1), n_cells=1
    )

    # Set the estimator values
    estimator_state.photo_ion_estimator[:] = photo_ion_array
    estimator_state.stim_recomb_estimator[:] = stim_recomb_array

    estimators_continuum = {
        "photoionization_rate_estimator": estimator_state.photo_ion_estimator,
        "stimulated_recombination_rate_estimator": estimator_state.stim_recomb_estimator,
    }

    (
        actual_photoionization_rate_coeff,
        actual_stimulated_recombination_rate_coeff,
    ) = solver.solve(
        estimators_continuum,
    )

    assert isinstance(actual_photoionization_rate_coeff, pd.DataFrame)
    assert actual_photoionization_rate_coeff.shape[0] == len(
        level2continuum_edge_idx
    )
    expected_photoionization_rate_coeff = pd.DataFrame(
        # These fixed estimator coefficients distinguish the raw 1e-5 channel
        # from the corresponding stimulated-recombination 1e-6 channel.
        [[0.0], [2.0e-5], [3.0e-5]],
        index=photoionization_index,
        columns=pd.Index([0], name="Shell No."),
    )
    pdt.assert_frame_equal(
        actual_photoionization_rate_coeff,
        expected_photoionization_rate_coeff,
        atol=0,
        rtol=1e-15,
    )
    expected_stimulated_recombination_rate_coeff = pd.DataFrame(
        [[0.0], [2.0e-6], [3.0e-6]],
        index=photoionization_index,
        columns=pd.Index([0], name="Shell No."),
    )
    pdt.assert_frame_equal(
        actual_stimulated_recombination_rate_coeff,
        expected_stimulated_recombination_rate_coeff,
        atol=0,
        rtol=1e-15,
    )


def test_estimated_solver_exposes_raw_coefficients_without_population_correction():
    photoionization_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    level2continuum_edge_idx = pd.Series([0], index=photoionization_index)
    solver = EstimatedPhotoionizationCoeffSolver(level2continuum_edge_idx)

    estimators_continuum = {
        # Unequal tallies make an unintended population correction observable.
        "photoionization_rate_estimator": pd.DataFrame(
            [[5.0]], index=photoionization_index, columns=[0]
        ),
        "stimulated_recombination_rate_estimator": pd.DataFrame(
            [[2.0]], index=photoionization_index, columns=[0]
        ),
    }
    photoionization, stimulated_recombination = solver.solve(
        estimators_continuum,
    )

    # Estimator-backed coefficients are raw inputs to rate assembly; applying
    # any population correction here would change these values.
    pdt.assert_frame_equal(
        photoionization,
        estimators_continuum["photoionization_rate_estimator"].set_axis(
            pd.Index([0], name="Shell No."), axis=1
        ),
    )
    pdt.assert_frame_equal(
        stimulated_recombination,
        estimators_continuum[
            "stimulated_recombination_rate_estimator"
        ].set_axis(pd.Index([0], name="Shell No."), axis=1),
    )


def test_photoionization_solver_does_not_add_hydrogen_policy_rows_to_helium(
    mock_photoionization_cross_sections,
):
    helium_index = pd.MultiIndex.from_tuples(
        [(2, 0, 0)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    helium_cross_sections = mock_photoionization_cross_sections.copy()
    helium_cross_sections.index = helium_index.repeat(2)
    solver = AnalyticPhotoionizationCoeffSolver(helium_cross_sections)
    radiation_field = DilutePlanckianRadiationField(
        np.array([10000.0, 10000.0]) * u.K, np.ones(2)
    )

    gamma, stimulated_recombination = solver.solve(
        radiation_field, np.array([10000.0, 10000.0]) * u.K
    )

    assert list(gamma.index) == [(2, 0, 0)]
    assert list(stimulated_recombination.index) == [(2, 0, 0)]


def test_lucy_recombination_factor_is_applied_once():
    index = pd.MultiIndex.from_tuples(
        [(1, 0, 1)],
        names=["atomic_number", "ion_number", "level_number"],
    )
    raw_coefficient = pd.DataFrame([[2.0]], index=index, columns=[0])
    phi_ik = pd.DataFrame([[3.0]], index=index, columns=[0])

    actual = apply_lte_level_to_ion_factor(raw_coefficient, phi_ik)

    # Recombination assembly owns this Lucy factor, and the raw coefficient
    # must remain reusable for the other shells and transitions.
    pdt.assert_frame_equal(
        actual,
        pd.DataFrame([[6.0]], index=index, columns=[0]),
    )
    pdt.assert_frame_equal(
        raw_coefficient,
        pd.DataFrame([[2.0]], index=index, columns=[0]),
    )
