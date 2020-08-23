import pytest
import numpy as np
from tardis.widgets.line_info import LineInfoWidget
from tardis.util.base import species_string_to_tuple


@pytest.fixture(scope="class")
def line_info_widget(simulation_verysimple):
    line_info_widget = LineInfoWidget.from_simulation(simulation_verysimple)
    # To attach event listeners to component widgets of line_info_widget
    _ = line_info_widget.display()
    return line_info_widget


@pytest.mark.parametrize(
    ("wavelength_range", "filter_mode"),
    [
        ([3000, 4000], "packet_out_nu"),
        ([3000, 4000], "packet_in_nu"),
        ([16200, 16300], "packet_out_nu"),
        (None, "packet_in_nu"),
    ],
)
class TestLineInfoWidgetData:
    """Tests for methods that handles data in LineInfoWidget."""

    def test_get_species_interactions(
        self, line_info_widget, wavelength_range, filter_mode
    ):
        species_interactions_df = line_info_widget.get_species_interactions(
            wavelength_range, filter_mode
        )

        if wavelength_range is None or wavelength_range == [16200, 16300]:
            # Dataframe contains all falsy values (proxy for empty)
            assert species_interactions_df.all(axis=None) == False
        else:
            # All values sum up to 1
            assert np.isclose(species_interactions_df.iloc[:, 0].sum(), 1)

            # Test shape of the dataframe
            expected_df_length = (
                line_info_widget.line_interaction_analysis[filter_mode]
                .last_line_in.groupby(["atomic_number", "ion_number"])
                .ngroups
            )
            assert species_interactions_df.shape == (expected_df_length, 1)

    @pytest.fixture
    def allowed_species(self, line_info_widget, wavelength_range, filter_mode):
        # Find species present within the selected wavelength range
        species_interactions_df = line_info_widget.get_species_interactions(
            wavelength_range, filter_mode
        )
        if species_interactions_df.all(axis=None) == False:
            allowed_species = None  # no species can be selected
        else:
            allowed_species = species_interactions_df.index
        return allowed_species

    @pytest.mark.parametrize("group_mode", ["exc", "de-exc", "both"])
    def test_get_last_line_counts(
        self, line_info_widget, allowed_species, filter_mode, group_mode
    ):
        """
        Test for get_last_line_counts() method.

        Since this method depends on get_species_interactions() so we need to
        make sure that we select only allowed species i.e. present within the
        wavelength range selected by the get_species_interactions()
        """
        if allowed_species is None:
            last_line_counts_df = line_info_widget.get_last_line_counts(
                None, filter_mode, group_mode
            )

            # Dataframe contains all falsy values (proxy for empty)
            assert last_line_counts_df.all(axis=None) == False

            return

        for selected_species in allowed_species:
            last_line_counts_df = line_info_widget.get_last_line_counts(
                selected_species, filter_mode, group_mode
            )

            last_lines_in = (
                line_info_widget.line_interaction_analysis[filter_mode]
                .last_line_in.xs(
                    key=species_string_to_tuple(selected_species),
                    level=["atomic_number", "ion_number"],
                )
                .reset_index()
            )

            last_lines_out = (
                line_info_widget.line_interaction_analysis[filter_mode]
                .last_line_out.xs(
                    key=species_string_to_tuple(selected_species),
                    level=["atomic_number", "ion_number"],
                )
                .reset_index()
            )

            if group_mode == "exc":
                expected_df_length = last_lines_in.groupby("line_id").ngroups
            elif group_mode == "de-exc":
                expected_df_length = last_lines_out.groupby("line_id").ngroups
            elif group_mode == "both":
                expected_df_length = last_lines_in.groupby(
                    ["line_id", last_lines_out["line_id"]]
                ).ngroups

            # Test shape of the dataframe
            assert last_line_counts_df.shape == (expected_df_length, 1)

