import pandas.testing as pdt
import pytest

from tardis.opacities.macro_atom.macroatom_solver import LegacyMacroAtomSolver
from tardis.opacities.macro_atom.macroatom_state import LegacyMacroAtomState
from tardis.opacities.opacity_solver import OpacitySolver


@pytest.fixture
def legacy_macro_atom_state(nb_simulation_verysimple):
    """Create a LegacyMacroAtomState for testing."""
    legacy_plasma = nb_simulation_verysimple.plasma

    opacity_solver = OpacitySolver(
        line_interaction_type="macroatom",
        disable_line_scattering=False,
    )
    opacity_state = opacity_solver.legacy_solve(legacy_plasma)

    macro_atom_state = LegacyMacroAtomSolver().solve(
        legacy_plasma.j_blues,
        legacy_plasma.atomic_data,
        opacity_state.tau_sobolev,
        legacy_plasma.stimulated_emission_factor,
        beta_sobolev=opacity_state.beta_sobolev,
    )
    # current implementation has 20 columns for the simulation

    return macro_atom_state


class TestLegacyMacroAtomStateSlicing:
    """Test class for LegacyMacroAtomState slicing functionality."""

    def test_getitem_slice(self, legacy_macro_atom_state):
        """Test __getitem__ with slice object."""
        original_state = legacy_macro_atom_state

        # Test slice from start
        sliced_state = original_state[0:3]

        assert isinstance(sliced_state, LegacyMacroAtomState)
        assert sliced_state.transition_probabilities.shape[1] == 3
        assert (
            sliced_state.transition_probabilities.shape[0]
            == original_state.transition_probabilities.shape[0]
        )

        # Verify sliced data matches original
        expected_columns = original_state.transition_probabilities.iloc[:, 0:3]
        pdt.assert_frame_equal(
            sliced_state.transition_probabilities, expected_columns
        )

        # Verify other attributes are copied correctly
        assert sliced_state.transition_type.equals(original_state.transition_type)

    def test_getitem_slice_middle(self, legacy_macro_atom_state):
        """Test __getitem__ with slice object from middle."""
        original_state = legacy_macro_atom_state

        # Test slice from middle (columns 5-10)
        sliced_state = original_state[5:10]

        assert isinstance(sliced_state, LegacyMacroAtomState)
        assert sliced_state.transition_probabilities.shape[1] == 5

        # Verify sliced data matches original
        expected_columns = original_state.transition_probabilities.iloc[:, 5:10]
        pdt.assert_frame_equal(
            sliced_state.transition_probabilities, expected_columns
        )

    def test_getitem_slice_with_step(self, legacy_macro_atom_state):
        """Test __getitem__ with slice object with step."""
        original_state = legacy_macro_atom_state

        # Test slice with step (every 2nd column)
        sliced_state = original_state[::2]

        assert isinstance(sliced_state, LegacyMacroAtomState)
        assert sliced_state.transition_probabilities.shape[1] == 10  # 20/2

        # Verify sliced data matches original
        expected_columns = original_state.transition_probabilities.iloc[:, ::2]
        pdt.assert_frame_equal(
            sliced_state.transition_probabilities, expected_columns
        )

    def test_getitem_single_integer(self, legacy_macro_atom_state):
        """Test __getitem__ with single integer index."""
        original_state = legacy_macro_atom_state

        # Test single column selection
        sliced_state = original_state[5]

        assert isinstance(sliced_state, LegacyMacroAtomState)
        assert sliced_state.transition_probabilities.shape[1] == 1

        # Verify sliced data matches original single column
        expected_column = original_state.transition_probabilities.iloc[:, 5:6]
        pdt.assert_frame_equal(
            sliced_state.transition_probabilities, expected_column
        )

    def test_slice_columns_copy_behavior(self, legacy_macro_atom_state):
        """Test _slice_columns with copy=True vs copy=False."""
        original_state = legacy_macro_atom_state

        # Test with copy=True
        sliced_copy = original_state._slice_columns(slice(0, 3), copy=True)
        assert sliced_copy.transition_type is not original_state.transition_type

        # Test with copy=False
        sliced_ref = original_state._slice_columns(slice(0, 3), copy=False)
        assert sliced_ref.transition_type is original_state.transition_type

    def test_edge_cases(self, legacy_macro_atom_state):
        """Test edge cases for slicing."""
        original_state = legacy_macro_atom_state

        # Test last column
        sliced_state = original_state[19]  # Last column (0-indexed)
        assert sliced_state.transition_probabilities.shape[1] == 1

        # Test empty slice
        sliced_state = original_state[5:5]  # Empty slice
        assert sliced_state.transition_probabilities.shape[1] == 0
        assert isinstance(sliced_state, LegacyMacroAtomState)

    def test_duck_typing_behavior(self, legacy_macro_atom_state):
        """Test that duck typing works correctly for different key types."""
        original_state = legacy_macro_atom_state

        # Test with slice object and integer
        sliced_with_slice = original_state[slice(0, 3)]
        sliced_with_int = original_state[0]

        assert isinstance(sliced_with_slice, LegacyMacroAtomState)
        assert isinstance(sliced_with_int, LegacyMacroAtomState)
        assert sliced_with_slice.transition_probabilities.shape[1] == 3
        assert sliced_with_int.transition_probabilities.shape[1] == 1

    @pytest.mark.parametrize("copy_flag", [True, False])
    def test_parametrized_copy_behavior(
        self, legacy_macro_atom_state, copy_flag
    ):
        """Parametrized test for copy behavior."""
        original_state = legacy_macro_atom_state
        sliced_state = original_state._slice_columns(
            slice(0, 5), copy=copy_flag
        )

        assert isinstance(sliced_state, LegacyMacroAtomState)
        assert sliced_state.transition_probabilities.shape[1] == 5

        # Verify data correctness regardless of copy mode
        expected_columns = original_state.transition_probabilities.iloc[:, 0:5]
        pdt.assert_frame_equal(
            sliced_state.transition_probabilities, expected_columns
        )

    def test_multiple_slicing_operations(self, legacy_macro_atom_state):
        """Test multiple consecutive slicing operations."""
        original_state = legacy_macro_atom_state

        # First slice (columns 0-9)
        first_slice = original_state[0:10]
        assert first_slice.transition_probabilities.shape[1] == 10

        # Second slice on the result (columns 2-5 of the first slice)
        second_slice = first_slice[2:6]
        assert second_slice.transition_probabilities.shape[1] == 4

        # Verify final result matches direct slicing
        direct_slice = original_state.transition_probabilities.iloc[:, 2:6]
        pdt.assert_frame_equal(
            second_slice.transition_probabilities, direct_slice
        )
