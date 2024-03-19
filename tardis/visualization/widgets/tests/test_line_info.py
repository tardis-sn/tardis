import pytest
import pandas as pd
import numpy as np
from plotly.callbacks import Points, BoxSelector
from tardis.visualization.widgets.line_info import LineInfoWidget
from tardis.util.base import species_string_to_tuple
from tardis.tests.test_util import monkeysession


@pytest.fixture(scope="class")
def line_info_widget(simulation_verysimple):
    line_info_widget = LineInfoWidget.from_simulation(simulation_verysimple)
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
        """
        Test for get_species_interactions() method.

        Checks shape of dataframe and whether all values sum up to 1 in cases
        where dataframe resulting dataframe should not be empty.
        """
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
        """
        For different combinations of wavelength_range and filter_mode
        parameters, it calls get_species_interactions on line_info_widget
        """
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
        wavelength range selected by the get_species_interactions(), which is
        being done here by allowed_species fixture.
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


class TestLineInfoWidgetEvents:
    """
    Test changes in table widgets data by triggering all possible events.

    This will make sure that all event listeners are working properly and
    updating data in tables accurately. The following four methods are to
    trigger each event (interaction) which is possible in LineInfoWidget.
    """

    @pytest.fixture(
        scope="class",
        params=[
            [2500, 3500],  # Wavelength range with plenty of line interactions
            [16200, 16300],  # Wavelength range with no line interactions
            None,  # No selection of wavelength range
        ],
    )
    def liw_with_selection(self, simulation_verysimple, request, monkeysession):
        """
        Makes different wavelength range selection on figure (specified by
        params) after creating a LineInfoWidget object.
        """
        liw = LineInfoWidget.from_simulation(simulation_verysimple)
        monkeysession.setattr(
            "tardis.visualization.widgets.line_info.is_notebook", lambda: True
        )
        # To attach event listeners to component widgets of line_info_widget
        _ = liw.display()

        selection_range = request.param

        # Since we cannot programatically make a Box selection on spectrum
        # so we have to directly call its event listener by passing
        # selected wavelength range in a BoxSelector object
        if selection_range:
            liw._spectrum_selection_handler(
                trace=liw.figure_widget.data[0],
                points=Points(),
                selector=BoxSelector(
                    xrange=selection_range,
                    yrange=[
                        -1.8e39,
                        1.8e40,
                    ],  # Not very relevant, approx height of box
                ),
            )

        return liw, selection_range

    def test_selection_on_plot(self, liw_with_selection):
        """
        Test if selection on spectrum plot, updates correct data in both
        the tables and total packets label.
        """
        # Since wavelength range selection is already made by liw_with_selection
        # fixture, we don't need to trigger selection event here again

        line_info_widget, selected_wavelength_range = liw_with_selection

        expected_species_interactions = (
            line_info_widget.get_species_interactions(
                wavelength_range=selected_wavelength_range,
                filter_mode=line_info_widget.FILTER_MODES[
                    line_info_widget.filter_mode_buttons.index
                ],
            )
        )

        pd.testing.assert_frame_equal(
            expected_species_interactions,
            line_info_widget.species_interactions_table.df,
        )

        expected_last_line_counts = line_info_widget.get_last_line_counts(
            selected_species=expected_species_interactions.index[0],
            filter_mode=line_info_widget.FILTER_MODES[
                line_info_widget.filter_mode_buttons.index
            ],
            group_mode=line_info_widget.GROUP_MODES[
                line_info_widget.group_mode_dropdown.index
            ],
        )

        pd.testing.assert_frame_equal(
            expected_last_line_counts,
            line_info_widget.last_line_counts_table.df,
        )

        if selected_wavelength_range in [None, [16200, 16300]]:
            expected_total_packets = 0
        else:
            expected_total_packets = expected_last_line_counts.iloc[:, 0].sum()
        assert expected_total_packets == int(
            line_info_widget.total_packets_label.widget.children[1].value
        )

    @pytest.mark.parametrize("selected_filter_mode_idx", [0, 1])
    def test_filter_mode_toggle(
        self,
        liw_with_selection,
        selected_filter_mode_idx,
    ):
        """
        Test if toggling filter_mode_buttons updates correct data in both
        the tables and total packets label.
        """
        line_info_widget, selected_wavelength_range = liw_with_selection

        # Toggle the filter_mode_buttons
        line_info_widget.filter_mode_buttons.index = selected_filter_mode_idx

        expected_species_interactions = (
            line_info_widget.get_species_interactions(
                wavelength_range=selected_wavelength_range,
                filter_mode=line_info_widget.FILTER_MODES[
                    selected_filter_mode_idx
                ],
            )
        )

        pd.testing.assert_frame_equal(
            expected_species_interactions,
            line_info_widget.species_interactions_table.df,
        )

        expected_last_line_counts = line_info_widget.get_last_line_counts(
            selected_species=expected_species_interactions.index[0],
            filter_mode=line_info_widget.FILTER_MODES[selected_filter_mode_idx],
            group_mode=line_info_widget.GROUP_MODES[
                line_info_widget.group_mode_dropdown.index
            ],
        )

        pd.testing.assert_frame_equal(
            expected_last_line_counts,
            line_info_widget.last_line_counts_table.df,
        )

        if selected_wavelength_range in [None, [16200, 16300]]:
            expected_total_packets = 0
        else:
            expected_total_packets = expected_last_line_counts.iloc[:, 0].sum()
        assert expected_total_packets == int(
            line_info_widget.total_packets_label.widget.children[1].value
        )

    def test_selection_on_species_intrctn_table(self, liw_with_selection):
        """
        Test if selection on each row in species_interaction_table updates
        correct data in last_line_counts_table and total packets label.
        """
        line_info_widget, _ = liw_with_selection

        for (
            selected_species
        ) in line_info_widget.species_interactions_table.df.index:
            # Select row in species_interactions_table
            line_info_widget.species_interactions_table.change_selection(
                [selected_species]
            )

            if bool(selected_species) == False:
                # When selected_species is a falsy value due to empty
                # species_interactions_table, use it as None in get_last_line_counts()
                selected_species = None

            expected_last_line_counts = line_info_widget.get_last_line_counts(
                selected_species=selected_species,
                filter_mode=line_info_widget.FILTER_MODES[
                    line_info_widget.filter_mode_buttons.index
                ],
                group_mode=line_info_widget.GROUP_MODES[
                    line_info_widget.group_mode_dropdown.index
                ],
            )

            pd.testing.assert_frame_equal(
                expected_last_line_counts,
                line_info_widget.last_line_counts_table.df,
            )

            if selected_species is None:
                expected_total_packets = 0
            else:
                expected_total_packets = expected_last_line_counts.iloc[
                    :, 0
                ].sum()
            assert expected_total_packets == int(
                line_info_widget.total_packets_label.widget.children[1].value
            )

    @pytest.mark.parametrize("selected_group_mode_idx", [0, 1, 2])
    def test_group_mode_change(
        self, liw_with_selection, selected_group_mode_idx
    ):
        """
        Test if selecting an option from group_mode_dropdown updates
        correct data in last_line_counts_table and total packets label.
        """
        line_info_widget, _ = liw_with_selection

        # Select the option in group_mode_dropdown
        line_info_widget.group_mode_dropdown.index = selected_group_mode_idx

        # For testing changes in last_line_counts_table data,
        # we're only considering the 1st row (0th index species)
        # in species_interactions_table
        if line_info_widget.last_line_counts_table.df.all(axis=None) == False:
            species0 = None
        else:
            species0 = line_info_widget.species_interactions_table.df.index[0]
            # Select 1st row in species_interaction_table, if not selected
            line_info_widget.species_interactions_table.change_selection(
                [species0]
            )

        expected_last_line_counts = line_info_widget.get_last_line_counts(
            selected_species=species0,
            filter_mode=line_info_widget.FILTER_MODES[
                line_info_widget.filter_mode_buttons.index
            ],
            group_mode=line_info_widget.GROUP_MODES[selected_group_mode_idx],
        )

        pd.testing.assert_frame_equal(
            expected_last_line_counts,
            line_info_widget.last_line_counts_table.df,
        )

        if species0 is None:
            expected_total_packets = 0
        else:
            expected_total_packets = expected_last_line_counts.iloc[:, 0].sum()
        assert expected_total_packets == int(
            line_info_widget.total_packets_label.widget.children[1].value
        )
