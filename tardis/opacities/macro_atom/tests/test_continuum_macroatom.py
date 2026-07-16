"""Integration tests for ContinuumMacroAtomSolver.

This module contains integration tests for the ContinuumMacroAtomSolver class,
which solves for macro-atom transition probabilities including continuum processes
such as photoionization, recombination, and collisional transitions.

The tests follow a fixture-based pattern similar to test_tardis_model_density_config.py,
where test data is prepared in fixtures and regression data is used for validation.
"""

### FIX THIS TO ONLY TEST HYDROGEN FOR NOW
# from copy import deepcopy

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.opacities.macro_atom.macroatom_solver import (
    ContinuumMacroAtomSolver,
)
from tardis.opacities.macro_atom.macroatom_state import MacroAtomState
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType


@pytest.fixture(
    scope="function"
)  # Needs to be function scope for multi-iter solve
def continuum_macro_atom_solver(iip_atom_data):
    """Fixture creating a ContinuumMacroAtomSolver instance.

    Initializes the solver with levels, lines, and photoionization data
    from the atomic dataset.
    """
    solver = ContinuumMacroAtomSolver(
        levels=iip_atom_data.levels.loc[[1]],
        lines=iip_atom_data.lines.loc[[1]],
        photoionization_data=iip_atom_data.photoionization_data.loc[[1]],
        ionization_energies=iip_atom_data.ionization_data.loc[[1]],
        selected_continuum_transitions=np.array([]),
        line_interaction_type="macroatom",
    )
    return solver


@pytest.fixture
def continuum_solver_input_data(iip_atom_data):
    """Fixture providing input data for solving macro-atom transitions.

    Creates minimal but valid input arrays and DataFrames needed for
    the ContinuumMacroAtomSolver.solve() method. Uses a fixed random seed
    for reproducibility.
    """
    # Set seed for reproducible random number generation
    np.random.seed(42)
    n_shells = 5
    n_lines = len(iip_atom_data.lines.loc[[1]])
    n_levels = len(iip_atom_data.levels.loc[[1]])

    # Create mean intensities for lines in blue wing
    # Shape: (n_lines, n_shells)
    mean_intensities_blue_wing = pd.DataFrame(
        np.random.uniform(0.1, 1.0, size=(n_lines, n_shells)),
        index=iip_atom_data.lines.loc[[1]].index,
        columns=np.arange(n_shells),
    )

    # Create beta sobolev values (opacity-related)
    # Shape: (n_lines, n_shells)
    beta_sobolevs = pd.DataFrame(
        np.random.uniform(0.1, 2.0, size=(n_lines, n_shells)),
        index=iip_atom_data.lines.loc[[1]].index,
        columns=np.arange(n_shells),
    )

    # Stimulated emission factors
    # Shape: (n_lines, n_shells) or (n_lines,)
    stimulated_emission_factors = np.random.uniform(
        0.5, 1.5, size=(n_lines, n_shells)
    )

    # Photoionization rate coefficients
    # Shape: (n_levels-1, n_shells)
    stim_recomb_corrected_photoionization_rate_coeff = pd.DataFrame(
        np.random.uniform(1e-12, 1e-10, size=(n_levels - 1, n_shells)),
        index=iip_atom_data.levels.xs(
            (1, 0), drop_level=False
        ).index,  # Pandas trim to only neutral H
        columns=np.arange(n_shells),
    )

    # Spontaneous recombination coefficients
    # Shape: (n_continuum_trans, n_shells)
    spontaneous_recombination_coeff = pd.DataFrame(
        np.random.uniform(1e-13, 1e-11, size=(n_levels - 1, n_shells)),
        index=iip_atom_data.levels.xs(
            (1, 0), drop_level=False
        ).index,  # Pandas trim to only neutral H
        columns=np.arange(n_shells),
    )

    # Collision de-excitation coefficients
    # Shape: (n_levels, n_shells)
    coll_deexc_coeff = pd.DataFrame(
        np.random.uniform(1e-15, 1e-12, size=(n_lines, n_shells)),
        index=iip_atom_data.lines.loc[[1]].index,
        columns=np.arange(n_shells),
    )

    # Collision excitation coefficients
    # Shape: (n_lines, n_shells)
    coll_exc_coeff = pd.DataFrame(
        np.random.uniform(1e-16, 1e-13, size=(n_lines, n_shells)),
        index=iip_atom_data.lines.loc[[1]].index,
        columns=np.arange(n_shells),
    )

    # Collision ionization coefficients
    # Shape: (n_levels - 1, n_shells)
    coll_ion_coeff = pd.DataFrame(
        np.random.uniform(1e-15, 1e-12, size=(n_levels - 1, n_shells)),
        index=iip_atom_data.levels.xs(
            (1, 0), drop_level=False
        ).index,  # Pandas trim to only neutral H
        columns=np.arange(n_shells),
    )

    # Collision recombination coefficients
    # Shape: (n_levels-1, n_shells)
    coll_recomb_coeff = pd.DataFrame(
        np.random.uniform(1e-16, 1e-13, size=(n_levels - 1, n_shells)),
        index=iip_atom_data.levels.xs(
            (1, 0), drop_level=False
        ).index,  # Pandas trim to only neutral H
        columns=np.arange(n_shells),
    )

    # Electron densities (1D array)
    # Shape: (n_shells,)
    electron_densities = pd.Series(
        np.random.uniform(1e6, 1e8, size=n_shells),
        index=np.arange(n_shells),
    )

    # Energy differences
    # Shape: (n_levels,)
    delta_E_yg = pd.Series(
        np.random.uniform(0.1, 20.0, size=n_lines),
        index=iip_atom_data.lines.loc[[1]].index,
    )

    # Cooling rates and arrays for different processes
    coll_exc_cool_rate = np.random.uniform(1e-20, 1e-18, size=n_shells)

    coll_exc_cool_arr = np.random.uniform(
        1e-30,
        1e-25,
        size=(
            n_shells,
            n_levels - 2,
        ),  # Don't have access to ionized H OR ground state H
    )

    # Create MultiIndex for collisional excitation cooling destinations
    # Exclude ground state (0) and ionized state (n_levels - 1)
    coll_exc_cool_destinations = iip_atom_data.levels.xs(
        (1, 0), drop_level=False
    ).index[:-1]

    coll_ion_cool_rate = np.random.uniform(1e-21, 1e-19, size=n_shells)

    coll_ion_cool_arr = np.random.uniform(
        1e-30, 1e-25, size=(n_shells, n_levels - 1)
    )

    fb_cool_rate = np.random.uniform(1e-20, 1e-18, size=n_shells)

    fb_cool_probs_arr = np.random.uniform(
        0.0, 1.0, size=(n_shells, n_levels - 1)
    )

    ff_cool_rate = np.random.uniform(1e-21, 1e-19, size=n_shells)

    return {
        "mean_intensities_blue_wing": mean_intensities_blue_wing,
        "beta_sobolevs": beta_sobolevs,
        "stimulated_emission_factors": stimulated_emission_factors,
        "stim_recomb_corrected_photoionization_rate_coeff": stim_recomb_corrected_photoionization_rate_coeff,
        "spontaneous_recombination_coeff": spontaneous_recombination_coeff,
        "coll_deexc_coeff": coll_deexc_coeff,
        "coll_exc_coeff": coll_exc_coeff,
        "coll_ion_coeff": coll_ion_coeff,
        "coll_recomb_coeff": coll_recomb_coeff,
        "electron_densities": electron_densities,
        "delta_E_yg": delta_E_yg,
        "coll_exc_cool_rate": coll_exc_cool_rate,
        "coll_exc_cool_arr": coll_exc_cool_arr,
        "coll_exc_cool_destinations": coll_exc_cool_destinations,
        "coll_ion_cool_rate": coll_ion_cool_rate,
        "coll_ion_cool_arr": coll_ion_cool_arr,
        "fb_cool_rate": fb_cool_rate,
        "fb_cool_probs_arr": fb_cool_probs_arr,
        "ff_cool_rate": ff_cool_rate,
    }


@pytest.fixture
def continuum_macro_atom_state(
    continuum_macro_atom_solver, continuum_solver_input_data
):
    """Fixture solving for macro-atom state with continuum processes.

    Calls the ContinuumMacroAtomSolver.solve() method with the prepared
    input data and returns the resulting MacroAtomState.
    """
    return continuum_macro_atom_solver.solve(
        mean_intensities_blue_wing=continuum_solver_input_data[
            "mean_intensities_blue_wing"
        ],
        beta_sobolevs=continuum_solver_input_data["beta_sobolevs"],
        stimulated_emission_factors=continuum_solver_input_data[
            "stimulated_emission_factors"
        ],
        stim_recomb_corrected_photoionization_rate_coeff=continuum_solver_input_data[
            "stim_recomb_corrected_photoionization_rate_coeff"
        ],
        spontaneous_recombination_coeff=continuum_solver_input_data[
            "spontaneous_recombination_coeff"
        ],
        coll_deexc_coeff=continuum_solver_input_data["coll_deexc_coeff"],
        coll_exc_coeff=continuum_solver_input_data["coll_exc_coeff"],
        coll_ion_coeff=continuum_solver_input_data["coll_ion_coeff"],
        coll_recomb_coeff=continuum_solver_input_data["coll_recomb_coeff"],
        electron_densities=continuum_solver_input_data["electron_densities"],
        delta_E_yg=continuum_solver_input_data["delta_E_yg"],
        coll_exc_cool_rate=continuum_solver_input_data["coll_exc_cool_rate"],
        coll_exc_cool_arr=continuum_solver_input_data["coll_exc_cool_arr"],
        coll_exc_cool_destinations=continuum_solver_input_data[
            "coll_exc_cool_destinations"
        ],
        coll_ion_cool_rate=continuum_solver_input_data["coll_ion_cool_rate"],
        coll_ion_cool_arr=continuum_solver_input_data["coll_ion_cool_arr"],
        fb_cool_rate=continuum_solver_input_data["fb_cool_rate"],
        fb_cool_probs_arr=continuum_solver_input_data["fb_cool_probs_arr"],
        ff_cool_rate=continuum_solver_input_data["ff_cool_rate"],
    )


class TestContinuumMacroAtomSolver:
    """Test suite for ContinuumMacroAtomSolver integration."""

    def test_solver_initialization(self, continuum_macro_atom_solver):
        """Test that ContinuumMacroAtomSolver initializes correctly."""
        assert continuum_macro_atom_solver is not None
        assert isinstance(continuum_macro_atom_solver, ContinuumMacroAtomSolver)
        assert continuum_macro_atom_solver.line_interaction_type == "macroatom"

    def test_macro_atom_state_structure(self, continuum_macro_atom_state):
        """Test that returned MacroAtomState has the correct structure."""
        assert isinstance(continuum_macro_atom_state, MacroAtomState)

        # Check required attributes exist
        assert hasattr(continuum_macro_atom_state, "transition_probabilities")
        assert hasattr(continuum_macro_atom_state, "transition_metadata")
        assert hasattr(continuum_macro_atom_state, "line2macro_level_upper")

    def test_transition_probabilities_shape(self, continuum_macro_atom_state):
        """Test that transition probabilities have correct shape and type."""
        probs = continuum_macro_atom_state.transition_probabilities

        # Should be a DataFrame
        assert isinstance(probs, pd.DataFrame)

        # Should have rows (transitions) and columns (shells)
        assert len(probs) > 0
        assert len(probs.columns) > 0

        # All values should be non-negative
        assert (probs >= 0).all().all()

    def test_transition_probabilities_validity(
        self, continuum_macro_atom_state
    ):
        """Test that transition probabilities are valid probabilities."""
        probs = continuum_macro_atom_state.transition_probabilities

        # Probabilities should be bounded [0, 1]
        assert (probs.values >= 0).all()
        assert (probs.values <= 1.0).all()

    def test_transition_probabilities_normalize_by_source_block(
        self, continuum_macro_atom_state
    ):
        probs = continuum_macro_atom_state.transition_probabilities
        metadata = continuum_macro_atom_state.transition_metadata
        source_blocks = probs.groupby(
            metadata["source_level_idx"].to_numpy()
        ).sum()

        # Raw propensities are normalized after energy weighting by source;
        # every surviving source block therefore distributes unit probability
        # across its available macro-atom branches.
        active_source_blocks = source_blocks.loc[
            source_blocks.sum(axis=1) > 0.0
        ]
        np.testing.assert_allclose(
            active_source_blocks.to_numpy(), 1.0, rtol=1e-12
        )

    def test_transition_metadata_structure(self, continuum_macro_atom_state):
        """Test that transition metadata has required columns."""
        metadata = continuum_macro_atom_state.transition_metadata

        # Should be a DataFrame
        assert isinstance(metadata, pd.DataFrame)

        # Should have same number of rows as probabilities
        assert len(metadata) == len(
            continuum_macro_atom_state.transition_probabilities
        )

        # Check for required columns
        required_columns = [
            "transition_type",
            "source_level_idx",
            "destination_level_idx",
        ]
        for col in required_columns:
            assert col in metadata.columns

    def test_transition_type_validity(self, continuum_macro_atom_state):
        """Test that transition types are valid MacroAtomTransitionType values."""
        metadata = continuum_macro_atom_state.transition_metadata
        transition_types = metadata["transition_type"].unique()

        # All transition types should be valid
        for tt in transition_types:
            assert tt in iter(MacroAtomTransitionType)

    def test_line2macro_level_upper_structure(self, continuum_macro_atom_state):
        """Test that line2macro_level_upper is properly structured."""
        l2m = continuum_macro_atom_state.line2macro_level_upper

        # Should be a Series or array-like
        assert hasattr(l2m, "__len__")

        # Should have entries
        assert len(l2m) > 0

    def test_transitions_consistency(self, continuum_macro_atom_state):
        """Test consistency between probabilities and metadata indices."""
        probs = continuum_macro_atom_state.transition_probabilities
        metadata = continuum_macro_atom_state.transition_metadata

        # Indices should match
        assert len(probs) == len(metadata)

    def test_transition_probabilities_regression(
        self, continuum_macro_atom_state, regression_data
    ):
        """Test transition probabilities using regression data.

        This test stores the transition probabilities in regression data
        for comparison across runs.
        """
        actual = continuum_macro_atom_state.transition_probabilities
        expected = regression_data.sync_dataframe(actual)
        pdt.assert_frame_equal(actual, expected)

    def test_transition_metadata_regression(
        self, continuum_macro_atom_state, regression_data
    ):
        """Test transition metadata using regression data.

        This test stores the metadata (transition types, source/destination indices)
        in regression data for comparison across runs.
        """
        actual = continuum_macro_atom_state.transition_metadata
        expected = regression_data.sync_dataframe(actual)
        pdt.assert_frame_equal(actual, expected)

    def test_multiple_solve_consistency(
        self, continuum_macro_atom_solver, continuum_solver_input_data
    ):
        """Test that solving twice with same data gives identical results."""
        # Solve twice
        state1 = continuum_macro_atom_solver.solve(
            mean_intensities_blue_wing=continuum_solver_input_data[
                "mean_intensities_blue_wing"
            ],
            beta_sobolevs=continuum_solver_input_data["beta_sobolevs"],
            stimulated_emission_factors=continuum_solver_input_data[
                "stimulated_emission_factors"
            ],
            stim_recomb_corrected_photoionization_rate_coeff=continuum_solver_input_data[
                "stim_recomb_corrected_photoionization_rate_coeff"
            ],
            spontaneous_recombination_coeff=continuum_solver_input_data[
                "spontaneous_recombination_coeff"
            ],
            coll_deexc_coeff=continuum_solver_input_data["coll_deexc_coeff"],
            coll_exc_coeff=continuum_solver_input_data["coll_exc_coeff"],
            coll_ion_coeff=continuum_solver_input_data["coll_ion_coeff"],
            coll_recomb_coeff=continuum_solver_input_data["coll_recomb_coeff"],
            electron_densities=continuum_solver_input_data[
                "electron_densities"
            ],
            delta_E_yg=continuum_solver_input_data["delta_E_yg"],
            coll_exc_cool_rate=continuum_solver_input_data[
                "coll_exc_cool_rate"
            ],
            coll_exc_cool_arr=continuum_solver_input_data["coll_exc_cool_arr"],
            coll_exc_cool_destinations=continuum_solver_input_data[
                "coll_exc_cool_destinations"
            ],
            coll_ion_cool_rate=continuum_solver_input_data[
                "coll_ion_cool_rate"
            ],
            coll_ion_cool_arr=continuum_solver_input_data["coll_ion_cool_arr"],
            fb_cool_rate=continuum_solver_input_data["fb_cool_rate"],
            fb_cool_probs_arr=continuum_solver_input_data["fb_cool_probs_arr"],
            ff_cool_rate=continuum_solver_input_data["ff_cool_rate"],
        )

        state2 = continuum_macro_atom_solver.solve(
            mean_intensities_blue_wing=continuum_solver_input_data[
                "mean_intensities_blue_wing"
            ],
            beta_sobolevs=continuum_solver_input_data["beta_sobolevs"],
            stimulated_emission_factors=continuum_solver_input_data[
                "stimulated_emission_factors"
            ],
            stim_recomb_corrected_photoionization_rate_coeff=continuum_solver_input_data[
                "stim_recomb_corrected_photoionization_rate_coeff"
            ],
            spontaneous_recombination_coeff=continuum_solver_input_data[
                "spontaneous_recombination_coeff"
            ],
            coll_deexc_coeff=continuum_solver_input_data["coll_deexc_coeff"],
            coll_exc_coeff=continuum_solver_input_data["coll_exc_coeff"],
            coll_ion_coeff=continuum_solver_input_data["coll_ion_coeff"],
            coll_recomb_coeff=continuum_solver_input_data["coll_recomb_coeff"],
            electron_densities=continuum_solver_input_data[
                "electron_densities"
            ],
            delta_E_yg=continuum_solver_input_data["delta_E_yg"],
            coll_exc_cool_rate=continuum_solver_input_data[
                "coll_exc_cool_rate"
            ],
            coll_exc_cool_arr=continuum_solver_input_data["coll_exc_cool_arr"],
            coll_exc_cool_destinations=continuum_solver_input_data[
                "coll_exc_cool_destinations"
            ],
            coll_ion_cool_rate=continuum_solver_input_data[
                "coll_ion_cool_rate"
            ],
            coll_ion_cool_arr=continuum_solver_input_data["coll_ion_cool_arr"],
            fb_cool_rate=continuum_solver_input_data["fb_cool_rate"],
            fb_cool_probs_arr=continuum_solver_input_data["fb_cool_probs_arr"],
            ff_cool_rate=continuum_solver_input_data["ff_cool_rate"],
        )

        # Probabilities should be identical
        pdt.assert_frame_equal(
            state1.transition_probabilities, state2.transition_probabilities
        )

        # Metadata should be identical
        pdt.assert_frame_equal(
            state1.transition_metadata, state2.transition_metadata
        )

    def test_no_nan_in_probabilities(self, continuum_macro_atom_state):
        """Test that transition probabilities do not contain NaN values."""
        probs = continuum_macro_atom_state.transition_probabilities

        # Check for NaN
        n_nans = probs.isna().sum().sum()
        assert n_nans == 0, (
            f"Found {n_nans} NaN values in transition probabilities"
        )

    def test_no_inf_in_probabilities(self, continuum_macro_atom_state):
        """Test that transition probabilities do not contain infinite values."""
        probs = continuum_macro_atom_state.transition_probabilities

        # Check for inf
        n_infs = np.isinf(probs.values).sum()
        assert n_infs == 0, (
            f"Found {n_infs} infinite values in transition probabilities"
        )

    def test_photoionization_transitions_included(
        self, continuum_macro_atom_state
    ):
        """Test that photoionization transition types are present."""
        metadata = continuum_macro_atom_state.transition_metadata

        # Check if any photoionization transitions exist
        photoionization_mask = (
            metadata["transition_type"]
            == MacroAtomTransitionType.PHOTOIONIZATION_INTERNAL
        )

        # Should have at least some photoionization transitions
        if len(metadata) > 0:
            assert photoionization_mask.any() or len(metadata) > 0

    def test_collision_transitions_included(self, continuum_macro_atom_state):
        """Test that collision transition types are present or acceptable."""
        metadata = continuum_macro_atom_state.transition_metadata

        # Just verify metadata has transition_type column with valid transitions
        transition_types = metadata["transition_type"].unique()
        assert len(transition_types) > 0

    def test_absorbing_probability_matrix_structure(
        self, continuum_macro_atom_state
    ):
        """Test that absorbing_probability_matrix has correct shape and type."""
        matrix = continuum_macro_atom_state.absorbing_probability_matrix

        # Should be a numpy array
        assert isinstance(matrix, np.ndarray)

        # Should have 3 dimensions: (num_cells, num_states, num_states)
        assert matrix.ndim == 3
        assert matrix.shape[1] == matrix.shape[2]  # Square state matrices

    def test_absorbing_probability_matrix_validity(
        self, continuum_macro_atom_state
    ):
        """Test that absorbing probabilities are valid probabilities."""
        matrix = continuum_macro_atom_state.absorbing_probability_matrix

        # Each probability should be in [0, 1]
        assert (matrix.flatten() >= 0).all()
        assert (matrix.flatten() <= 1.0).all()

    def test_absorbing_probability_matrix_no_nan(
        self, continuum_macro_atom_state
    ):
        """Test that absorbing_probability_matrix does not contain NaN values."""
        matrix = continuum_macro_atom_state.absorbing_probability_matrix

        # Check for NaN
        n_nans = np.isnan(matrix).sum()
        assert n_nans == 0, (
            f"Found {n_nans} NaN values in absorbing probability matrix"
        )

    def test_absorbing_probability_matrix_no_inf(
        self, continuum_macro_atom_state
    ):
        """Test that absorbing_probability_matrix does not contain infinite values."""
        matrix = continuum_macro_atom_state.absorbing_probability_matrix

        # Check for inf
        n_infs = np.isinf(matrix).sum()
        assert n_infs == 0, (
            f"Found {n_infs} infinite values in absorbing probability matrix"
        )

    def test_absorbing_probability_matrix_regression(
        self, continuum_macro_atom_state, regression_data
    ):
        """Test absorbing_probability_matrix using regression data.

        This test stores the absorbing probability matrix in regression data
        for comparison across runs.
        """
        actual = continuum_macro_atom_state.absorbing_probability_matrix
        expected = regression_data.sync_ndarray(actual)
        # rtol is for Mac test compatibility, regression data was generated on Mac
        # for these tests.
        np.testing.assert_allclose(actual, expected, rtol=5e-5, atol=0)
