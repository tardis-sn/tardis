"""Tests for the Grotrian Widget"""

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from tardis.visualization.widgets.grotrian import (
    GrotrianPlot,
    GrotrianWidget,
)


@pytest.fixture(scope="module")
def grotrian_plot(simulation_simple_tracked):
    """Creates a GrotrianPlot from the tracked simulation fixture.

    Parameters
    ----------
    simulation_simple_tracked : tardis.simulation.base.Simulation
        Simulation object.

    Returns
    -------
    tardis.visualization.widgets.grotrian.GrotrianPlot
    """
    grotrian_plot = GrotrianPlot.from_simulation(simulation_simple_tracked)

    line_analysis = grotrian_plot._line_interaction_analysis[
        grotrian_plot.filter_mode
    ]
    species_group = line_analysis.last_line_in.groupby(
        ["atomic_number", "ion_number"]
    )
    first_species = list(species_group.groups.keys())[0]
    atomic_number, ion_number = first_species

    grotrian_plot._atomic_number = atomic_number
    grotrian_plot._ion_number = ion_number


class TestGrotrianPlot:
    """Tests for Grotrian Plot Class."""

    def test_compute_transitions(self, grotrian_plot, regression_data):
        grotrian_plot._compute_transitions()

        excitation_lines = grotrian_plot.excite_lines
        deexcitation_lines = grotrian_plot.deexcite_lines

        expected_excitation_lines = regression_data.sync_dataframe(excitation_lines)
        expected_deexcitation_lines = regression_data.sync_dataframe(deexcitation_lines)

        pdt.assert_frame_equal(excitation_lines, expected_excitation_lines)
        pdt.assert_frame_equal(deexcitation_lines, expected_deexcitation_lines)

    def test_compute_level_data(self, grotrian_plot, regression_data):
        grotrian_plot._compute_level_data()

        level_data = grotrian_plot.level_data

        expected_level_data = regression_data.sync_dataframe(level_data)

        pdt.assert_frame_equal(level_data, expected_level_data)


@pytest.fixture(scope="module")
def grotrian_widget(simulation_simple_tracked, monkeysession):
    """Creates a GrotrianWidget from the tracked simulation fixture.

    Returns
    -------
    GrotrianWidget
        A GrotrianWidget instance with display() called.
    """
    widget = GrotrianWidget.from_simulation(simulation_simple_tracked)
    monkeysession.setattr(
        "tardis.util.environment.Environment.is_notebook", lambda: True
    )
    widget.display()
    return widget


class TestGrotrianWidgetEvents:
    """Tests for GrotrianWidget event handlers triggered by widget changes."""

    def test_ion_selector_updates_plot_species(self, grotrian_widget):
        """Changing the ion_selector dropdown should update the GrotrianPlot's
        atomic_number and ion_number to match the selected species."""
        species_options = grotrian_widget.ion_selector.options
        if len(species_options) < 2:
            pytest.skip("Not enough species to test ion change")

        # Switch to the second available species
        grotrian_widget.ion_selector.value = species_options[1]

        from tardis.util.base import species_string_to_tuple

        expected_atomic, expected_ion = species_string_to_tuple(
            species_options[1]
        )
        assert grotrian_widget.plot._atomic_number == expected_atomic
        assert grotrian_widget.plot._ion_number == expected_ion

    def test_shell_selector_updates_plot_shell(self, grotrian_widget):
        """Changing the shell_selector dropdown should update the GrotrianPlot's
        shell attribute."""
        # Select a specific shell
        grotrian_widget.shell_selector.value = "1"
        assert grotrian_widget.plot.shell == 1

    def test_shell_selector_all_sets_none(self, grotrian_widget):
        """Selecting 'All' in the shell_selector should set the GrotrianPlot's
        shell attribute to None (averaging across shells)."""
        grotrian_widget.shell_selector.value = "All"
        assert grotrian_widget.plot.shell is None

    def test_max_level_selector_updates_plot(self, grotrian_widget):
        """Changing the max_level_selector should update the GrotrianPlot's
        max_levels attribute and recompute level data."""
        original_max = grotrian_widget.max_level_selector.value
        grotrian_widget.max_level_selector.value = 5
        assert grotrian_widget.plot.max_levels == 5
        assert len(grotrian_widget.plot.level_data) <= 6
        # Restore
        grotrian_widget.max_level_selector.value = original_max

    def test_y_scale_toggle_updates_plot(self, grotrian_widget):
        """Toggling the y_scale_selector should update the GrotrianPlot's
        y_scale attribute."""
        grotrian_widget.y_scale_selector.value = "Linear"
        assert grotrian_widget.plot.y_scale == "Linear"

        grotrian_widget.y_scale_selector.value = "Log"
        assert grotrian_widget.plot.y_scale == "Log"

    def test_wavelength_slider_updates_plot(self, grotrian_widget):
        """Changing the wavelength_range_selector should update the GrotrianPlot's
        min_wavelength and max_wavelength."""
        min_wl = grotrian_widget.wavelength_range_selector.min
        max_wl = grotrian_widget.wavelength_range_selector.max

        # Set to a sub-range
        mid = (min_wl + max_wl) / 2
        grotrian_widget.wavelength_range_selector.value = (min_wl, mid)

        assert grotrian_widget.plot.min_wavelength == pytest.approx(
            min_wl, rel=0.1
        )
        # max_wavelength gets +1 in the handler
        assert grotrian_widget.plot.max_wavelength == pytest.approx(
            mid + 1, rel=0.1
        )

    def test_wavelength_resetter_on_ion_change(self, grotrian_widget):
        """When the ion is changed, the wavelength slider should be reset
        to the new wavelength range of the selected species."""
        species_options = grotrian_widget.ion_selector.options

        # Trigger the ion change
        grotrian_widget.ion_selector.value = species_options[0]

        # The slider min/max should match the plot's wavelength range
        if grotrian_widget.plot.min_wavelength is not None:
            assert (
                grotrian_widget.wavelength_range_selector.min
                == grotrian_widget.plot.min_wavelength
            )
            assert (
                grotrian_widget.wavelength_range_selector.max
                == grotrian_widget.plot.max_wavelength
            )
