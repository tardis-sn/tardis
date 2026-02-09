import pandas as pd
import pytest

from tardis.opacities.macro_atom.macroatom_continuum_transitions import (
    collisional_transition_deexc_to_k_packet,
    collisional_transition_excitation_cool,
    collisional_transition_internal_down,
    collisional_transition_excitation_to_macro_atom,
    continuum_adiabatic_cooling,
    continuum_free_free_cooling,
    continuum_transition_photoionization,
    continuum_transition_recombination_emission,
    continuum_transition_recombination_internal,
    probability_adiabatic_cooling,
    probability_collision_deexc_to_k_packet,
    probability_collision_excitation_cool,
    probability_collision_internal_down,
    probability_collision_exc_to_macro,
    probability_free_free_cooling,
    probability_photoionization,
    probability_recombination_emission,
    probability_recombination_internal,
)
from tardis.opacities.macro_atom.macroatom_line_transitions import CONST_H_CGS
from tardis.transport.montecarlo.macro_atom import MacroAtomTransitionType


class TestRecombinationTransitions:
    """Test recombination-related transition functions."""

    @pytest.fixture
    def sample_recomb_data(self):
        """Create sample data for recombination tests."""
        # Create MultiIndex for atomic levels (atomic_number, ion_number, level_number)
        index = pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (2, 0, 0), (2, 1, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        )

        # Sample recombination coefficients (temperature-dependent)
        recomb_coeff = pd.DataFrame(
            [[1e-12, 2e-12], [3e-12, 4e-12], [5e-12, 6e-12], [7e-12, 8e-12]],
            index=index,
            columns=[0, 1],  # shell indices
        )

        # Sample level energies
        level_energies = pd.Series(
            [0.0, 10.2, 0.0, 24.6],  # eV
            index=index,
        )

        # Sample frequencies
        frequencies = pd.Series(
            [2.18e15, 3.29e15, 5.45e15, 1.31e16],  # Hz
            index=index,
        )

        return recomb_coeff, level_energies, frequencies

    def test_probability_recombination_internal(self, sample_recomb_data):
        """Test internal recombination probability calculation."""
        recomb_coeff, level_energies, _ = sample_recomb_data

        result = probability_recombination_internal(
            recomb_coeff, level_energies
        )

        # Check shape and type
        assert isinstance(result, pd.DataFrame)
        assert result.shape == recomb_coeff.shape
        assert result.index.equals(recomb_coeff.index)

        # Check values - should be recomb_coeff * level_energies
        expected = recomb_coeff.multiply(level_energies, axis=0)
        pd.testing.assert_frame_equal(result, expected)

        # Check that ground states (energy=0) give zero probabilities
        ground_state_mask = level_energies == 0.0
        assert (result.loc[ground_state_mask] == 0.0).all().all()

    def test_continuum_transition_recombination_internal(
        self, sample_recomb_data
    ):
        """Test internal recombination transition with metadata."""
        recomb_coeff, level_energies, _ = sample_recomb_data

        probabilities, metadata = continuum_transition_recombination_internal(
            recomb_coeff, level_energies
        )

        # Test probabilities
        assert isinstance(probabilities, pd.DataFrame)
        assert "source" in probabilities.columns

        # Test metadata
        assert isinstance(metadata, pd.DataFrame)
        required_cols = [
            "transition_line_id",
            "source",
            "destination",
            "transition_type",
            "transition_line_idx",
            "photoionization_key_idx",
        ]
        for col in required_cols:
            assert col in metadata.columns

        # Check transition type
        assert (
            metadata.transition_type == MacroAtomTransitionType.RECOMB_INTERNAL
        ).all()

        # Check that metadata index matches probabilities
        assert metadata.index.equals(probabilities.index)

    def test_probability_recombination_emission(self, sample_recomb_data):
        """Test emission recombination probability calculation."""
        recomb_coeff, _, frequencies = sample_recomb_data

        result = probability_recombination_emission(recomb_coeff, frequencies)

        # Check shape and type
        assert isinstance(result, pd.DataFrame)
        assert result.shape == recomb_coeff.shape

        # Check values include Planck constant factor

        expected = recomb_coeff.multiply(frequencies, axis=0) * CONST_H_CGS
        pd.testing.assert_frame_equal(result, expected)

    def test_continuum_transition_recombination_emission(
        self, sample_recomb_data
    ):
        """Test emission recombination transition with metadata."""
        recomb_coeff, _, frequencies = sample_recomb_data

        probabilities, metadata = continuum_transition_recombination_emission(
            recomb_coeff, frequencies
        )

        # Test metadata has correct transition type
        assert (
            metadata.transition_type == MacroAtomTransitionType.BF_EMISSION
        ).all()

        # Test that source column is added to probabilities
        assert "source" in probabilities.columns


class TestPhotoionizationTransitions:
    """Test photoionization-related transition functions."""

    @pytest.fixture
    def sample_photoionization_data(self):
        """Create sample photoionization data."""
        index = pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (2, 0, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        )

        rate_coeff = pd.DataFrame(
            [[1e-18, 2e-18], [3e-18, 4e-18], [5e-18, 6e-18]],
            index=index,
            columns=[0, 1],
        )

        level_energies = pd.Series([13.6, 3.4, 24.6], index=index)

        return rate_coeff, level_energies

    def test_probability_photoionization(self, sample_photoionization_data):
        """Test photoionization probability calculation."""
        rate_coeff, level_energies = sample_photoionization_data

        result = probability_photoionization(rate_coeff, level_energies)

        # Check basic properties
        assert isinstance(result, pd.DataFrame)
        assert result.shape == rate_coeff.shape

        # Check calculation
        expected = rate_coeff.multiply(level_energies, axis=0)
        pd.testing.assert_frame_equal(result, expected)

    def test_continuum_transition_photoionization(
        self, sample_photoionization_data
    ):
        """Test photoionization transition with metadata."""
        rate_coeff, level_energies = sample_photoionization_data

        probabilities, metadata = continuum_transition_photoionization(
            rate_coeff, level_energies
        )

        # Check transition type
        assert (
            metadata.transition_type == MacroAtomTransitionType.PHOTOIONIZATION
        ).all()

        # Check metadata structure
        assert len(metadata) == len(probabilities)


class TestCoolingTransitions:
    """Test cooling-related transition functions."""

    @pytest.fixture
    def sample_cooling_data(self):
        """Create sample data for cooling tests."""
        # Electron densities and temperatures
        electron_densities = pd.Series([1e10, 2e10, 3e10], index=[0, 1, 2])
        t_electrons = pd.Series([5000, 6000, 7000], index=[0, 1, 2])
        time_explosion = 1e6  # seconds

        # Ion number density with MultiIndex (for free-free cooling)
        ion_index = pd.MultiIndex.from_tuples(
            [(1, 0), (1, 1), (2, 0), (2, 1)],
            names=["atomic_number", "ion_number"],
        )
        ion_densities = pd.DataFrame(
            [
                [1e8, 2e8, 3e8],
                [5e7, 1e8, 1.5e8],
                [1e7, 2e7, 3e7],
                [5e6, 1e7, 1.5e7],
            ],
            index=ion_index,
            columns=[0, 1, 2],
        )

        return electron_densities, t_electrons, time_explosion, ion_densities

    def test_probability_adiabatic_cooling(self, sample_cooling_data):
        """Test adiabatic cooling probability calculation (currently not implemented)."""
        electron_densities, t_electrons, time_explosion, _ = sample_cooling_data

        # Function currently raises NotImplementedError
        with pytest.raises(NotImplementedError):
            probability_adiabatic_cooling(
                electron_densities, t_electrons, time_explosion
            )

    def test_continuum_adiabatic_cooling(self, sample_cooling_data):
        """Test adiabatic cooling with metadata (currently not implemented)."""
        electron_densities, t_electrons, time_explosion, _ = sample_cooling_data

        # Function currently raises NotImplementedError due to probability function
        with pytest.raises(NotImplementedError):
            continuum_adiabatic_cooling(
                electron_densities, t_electrons, time_explosion
            )

    def test_probability_free_free_cooling(self, sample_cooling_data):
        """Test free-free cooling probability calculation (currently not implemented)."""
        electron_densities, t_electrons, _, ion_densities = sample_cooling_data

        # Function currently raises NotImplementedError
        with pytest.raises(NotImplementedError):
            probability_free_free_cooling(
                ion_densities, electron_densities, t_electrons
            )

    def test_continuum_free_free_cooling(self, sample_cooling_data):
        """Test free-free cooling with metadata (currently not implemented)."""
        electron_densities, t_electrons, _, ion_densities = sample_cooling_data

        # Function currently raises NotImplementedError due to probability function
        with pytest.raises(NotImplementedError):
            continuum_free_free_cooling(
                ion_densities, electron_densities, t_electrons
            )


class TestCollisionalTransitions:
    """Test collisional transition functions."""

    @pytest.fixture
    def sample_collisional_data(self):
        """Create sample data for collisional tests."""
        # Collisional coefficients with MultiIndex
        coll_index = pd.MultiIndex.from_tuples(
            [(1, 0, 0, 1), (1, 0, 1, 2), (2, 0, 0, 1)],
            names=[
                "atomic_number",
                "ion_number",
                "level_number_lower",
                "level_number_upper",
            ],
        )

        coll_coeff = pd.DataFrame(
            [[1e-8, 2e-8, 3e-8], [4e-8, 5e-8, 6e-8], [7e-8, 8e-8, 9e-8]],
            index=coll_index,
            columns=[0, 1, 2],
        )

        electron_densities = pd.Series([1e10, 2e10, 3e10], index=[0, 1, 2])

        # Energy differences
        delta_E = pd.Series([10.2, 12.1, 24.6], index=coll_index)

        # Lower level energies
        energy_lower_index = pd.MultiIndex.from_tuples(
            [(1, 0, 0), (1, 0, 1), (2, 0, 0)],
            names=["atomic_number", "ion_number", "level_number"],
        )
        energy_lowers = pd.Series([0.0, 10.2, 0.0], index=energy_lower_index)

        # Level number densities
        level_densities = pd.Series(
            [1e12, 5e11, 2e12], index=energy_lower_index
        )

        # Mock atom_data object
        class MockAtomData:
            def __init__(self):
                self.levels = pd.DataFrame({"energy": energy_lowers})

        atom_data = MockAtomData()

        return (
            coll_coeff,
            electron_densities,
            delta_E,
            energy_lowers,
            level_densities,
            atom_data,
        )

    def test_probability_collision_deexc_to_k_packet(
        self, sample_collisional_data
    ):
        """Test collisional de-excitation to k-packet probability."""
        coll_coeff, electron_densities, delta_E, _, _, _ = (
            sample_collisional_data
        )

        result = probability_collision_deexc_to_k_packet(
            coll_coeff, electron_densities, delta_E
        )

        # Check that result is computed correctly
        assert isinstance(result, pd.DataFrame)

        # Check that energies are properly weighted (function doesn't do groupby)
        expected = (coll_coeff * electron_densities).multiply(
            delta_E.values, axis=0
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_collisional_transition_deexc_to_k_packet(
        self, sample_collisional_data
    ):
        """Test collisional de-excitation transition with metadata."""
        coll_coeff, electron_densities, delta_E, _, _, _ = (
            sample_collisional_data
        )

        probabilities, metadata = collisional_transition_deexc_to_k_packet(
            coll_coeff, electron_densities, delta_E
        )

        # Check metadata
        assert (
            metadata.transition_type
            == MacroAtomTransitionType.COLL_DOWN_TO_K_PACKET
        ).all()
        # Destination is actually a tuple (k, -99, -99), check just the first element
        assert all(dest[0] == "k" for dest in metadata.destination)

    def test_probability_collision_internal_down(self, sample_collisional_data):
        """Test collisional internal down probability."""
        coll_coeff, electron_densities, _, energy_lowers, _, _ = (
            sample_collisional_data
        )

        result = probability_collision_internal_down(
            coll_coeff, electron_densities, energy_lowers
        )

        # Check calculation
        expected = (coll_coeff * electron_densities).multiply(
            energy_lowers.values, axis=0
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_collisional_transition_internal_down(
        self, sample_collisional_data
    ):
        """Test collisional internal down transition with metadata."""
        coll_coeff, electron_densities, _, energy_lowers, _, _ = (
            sample_collisional_data
        )

        probabilities, metadata = collisional_transition_internal_down(
            coll_coeff, electron_densities, energy_lowers
        )

        # Check metadata
        assert (
            metadata.transition_type
            == MacroAtomTransitionType.COLL_DOWN_INTERNAL
        ).all()

    def test_probability_collision_internal_up(self, sample_collisional_data):
        """Test collisional internal up probability."""
        coll_coeff, electron_densities, _, energy_lowers, _, _ = (
            sample_collisional_data
        )

        result = probability_collision_exc_to_macro(
            coll_coeff, electron_densities, energy_lowers
        )

        # Check calculation
        expected = (coll_coeff * electron_densities).multiply(
            energy_lowers.values, axis=0
        )
        pd.testing.assert_frame_equal(result, expected)

    def test_collisional_transition_internal_up(self, sample_collisional_data):
        """Test collisional internal up transition with metadata."""
        coll_coeff, electron_densities, _, energy_lowers, _, _ = (
            sample_collisional_data
        )

        probabilities, metadata = (
            collisional_transition_excitation_to_macro_atom(
                coll_coeff, electron_densities, energy_lowers
            )
        )

        # Check metadata
        assert (
            metadata.transition_type
            == MacroAtomTransitionType.COLL_EXC_COOL_TO_MACRO
        ).all()

    def test_probability_collision_excitation_cool(
        self, sample_collisional_data
    ):
        """Test collisional excitation cooling probability."""
        coll_coeff, electron_densities, delta_E, _, level_densities, _ = (
            sample_collisional_data
        )

        # Create compatible lower indices for the test
        # Use MultiIndex tuples matching the level_densities structure
        lower_indices = coll_coeff.index.droplevel("level_number_upper")

        result = probability_collision_excitation_cool(
            coll_coeff,
            electron_densities,
            delta_E,
            level_densities,
            lower_indices,
        )

        # Check that result is computed properly
        assert isinstance(result, pd.DataFrame)
        # Just check basic shape and that it computed without error
        assert result.shape[1] == coll_coeff.shape[1]  # Same number of columns

    def test_collisional_transition_excitation_cool(
        self, sample_collisional_data
    ):
        """Test collisional excitation cooling transition with metadata."""
        coll_coeff, electron_densities, delta_E, _, level_densities, _ = (
            sample_collisional_data
        )

        # Create compatible lower indices that match level_densities index
        lower_indices = (
            level_densities.index
        )  # Use the actual level_densities index

        probabilities, metadata = collisional_transition_excitation_cool(
            coll_coeff,
            electron_densities,
            delta_E,
            level_densities,
            lower_indices,
        )

        # Check metadata
        assert (
            metadata.transition_type
            == MacroAtomTransitionType.COLL_EXC_COOL_TO_MACRO
        ).all()
        # Source is actually a tuple (k, -99, -99), check just the first element
        assert all(src[0] == "k" for src in metadata.source)


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_dataframes(self):
        """Test behavior with empty DataFrames."""
        empty_df = pd.DataFrame()
        empty_series = pd.Series([], dtype=float)

        # Test that functions handle empty inputs gracefully
        result = probability_recombination_internal(empty_df, empty_series)
        assert result.empty

        result = probability_recombination_emission(empty_df, empty_series)
        assert result.empty

    def test_zero_values(self):
        """Test behavior with zero values."""
        # Create data with zeros
        index = pd.MultiIndex.from_tuples(
            [(1, 0, 0)], names=["atomic_number", "ion_number", "level_number"]
        )
        zero_df = pd.DataFrame([[0.0]], index=index, columns=[0])
        zero_series = pd.Series([0.0], index=index)

        result = probability_recombination_internal(zero_df, zero_series)
        assert (result == 0.0).all().all()

    def test_negative_values_handling(self):
        """Test that functions handle negative inputs appropriately."""
        index = pd.MultiIndex.from_tuples(
            [(1, 0, 0)], names=["atomic_number", "ion_number", "level_number"]
        )

        # Negative coefficients should be allowed (could represent different physics)
        neg_df = pd.DataFrame([[-1e-12]], index=index, columns=[0])
        pos_series = pd.Series([10.0], index=index)

        result = probability_recombination_internal(neg_df, pos_series)
        assert (result < 0).all().all()  # Result should be negative
