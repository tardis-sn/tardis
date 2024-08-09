import numpy as np
import pytest
from numpy import testing as npt
from copy import deepcopy

from tardis.base import run_tardis
from tardis.visualization.tools.liv_plot import LIVPlotter
from tardis.tests.fixtures.regression_data import RegressionData


def make_valid_name(testid):
    """
    Sanitize pytest IDs to make them valid HDF group names.

    Parameters
    ----------
    testid : str
        ID to sanitize.

    Returns
    -------
    testid : str
        Sanitized ID.
    """
    return "_" + testid.replace("-", "_")


def idfn(val):
    return "_".join(f"{k}={v}" for k, v in val.items())


def pad_list(lst, length, pad_value=None):
    """
    Pad a list to a specified length with a given value.

    Parameters
    ----------
    lst : list
        List to pad.
    length : int
        Desired length.
    pad_value : any
        Value to pad the list with.

    Returns
    -------
    list
        Padded list.
    """
    return lst + [pad_value] * (length - len(lst))


@pytest.fixture(scope="module")
def simulation_simple(config_verysimple, atomic_dataset):
    """
    Run a simple TARDIS simulation for testing.

    Parameters
    ----------
    config_verysimple : tardis.io.config_reader.Configuration
        Configuration object for a very simple simulation.
    atomic_dataset : str or tardis.atomic.AtomData
        Atomic data.

    Returns
    -------
    sim: tardis.simulation.base.Simulation
        Simulation object.
    """
    config_verysimple.montecarlo.iterations = 3
    config_verysimple.montecarlo.no_of_packets = 4000
    config_verysimple.montecarlo.last_no_of_packets = -1
    config_verysimple.spectrum.virtual.virtual_packet_logging = True
    config_verysimple.montecarlo.no_of_virtual_packets = 1
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config_verysimple,
        atom_data=atomic_data,
    )
    return sim


@pytest.fixture(scope="class")
def plotter(simulation_simple):
    """
    Create a LIVPlotter object.

    Parameters
    ----------
    simulation_simple : tardis.simulation.base.Simulation
        Simulation object.

    Returns
    -------
    tardis.visualization.tools.liv_plot.LIVPlotter
    """
    return LIVPlotter.from_simulation(simulation_simple)


class TestLIVPlotter:
    """Test the LIVPlotter class."""

    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"]]
    )
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("nelements", [1, None])
    def test_parse_species_list(
        self,
        plotter,
        species_list,
        packets_mode,
        nelements,
        regression_data,
    ):
        """
        Test _parse_species_list method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.liv_plot.LIVPlotter
        species_list : List of species to plot.
        packets_mode : str, Packet mode, either 'virtual' or 'real'.
        nelements : int, Number of elements to include in plot.
        """
        plotter._parse_species_list(
            species_list=species_list,
            packets_mode=packets_mode,
            nelements=nelements,
        )

        actual_species_list = plotter._species_list
        actual_keep_colour = plotter._keep_colour
        flat_list = [
            item
            for sublist in list(plotter._species_mapped.values())
            for item in sublist
        ]
        actual_species_mapped = flat_list
        actual = (
            actual_species_list + actual_keep_colour + actual_species_mapped
        )
        expected = regression_data.sync_ndarray(actual)
        npt.assert_array_equal(actual, expected)

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"]]
    )
    @pytest.mark.parametrize("cmapname", ["jet"])
    @pytest.mark.parametrize("num_bins", [10, 25])
    @pytest.mark.parametrize("nelements", [1, None])
    def test_prepare_plot_data(
        self,
        plotter,
        packets_mode,
        species_list,
        cmapname,
        num_bins,
        nelements,
        regression_data,
    ):
        """
        Test _prepare_plot_data method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.liv_plot.LIVPlotter
        packets_mode : str, Packet mode, either 'virtual' or 'real'.
        species_list : list of species to plot
        cmapname : str
        num_bins : int, Number of bins for regrouping within the same range.
        nelements : int, Number of elements to include in plot.
        """
        plotter._prepare_plot_data(
            packets_mode=packets_mode,
            species_list=species_list,
            cmapname=cmapname,
            num_bins=num_bins,
            nelements=nelements,
        )
        plot_data_numeric = [
            [q.value for q in row] for row in plotter.plot_data
        ]
        flat_list = [item for sublist in plot_data_numeric for item in sublist]

        actual_plot_data_list = flat_list
        flat_list = [
            item for sublist in plotter.plot_colors for item in sublist
        ]
        actual_plot_colors = flat_list
        actual_new_bin_edges = list(plotter.new_bin_edges.value)
        actual = np.array(
            actual_plot_data_list + actual_plot_colors + actual_new_bin_edges
        )
        expected = regression_data.sync_ndarray(actual)
        npt.assert_array_equal(actual, expected)

    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"], None]
    )
    @pytest.mark.parametrize("nelements", [1, None])
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("xlog_scale", [True, False])
    @pytest.mark.parametrize("ylog_scale", [True, False])
    @pytest.mark.parametrize("num_bins", [10, 25])
    @pytest.mark.parametrize("velocity_range", [(18000, 25000)])
    def test_generate_plot_mpl(
        self,
        plotter,
        species_list,
        nelements,
        packets_mode,
        xlog_scale,
        ylog_scale,
        num_bins,
        velocity_range,
        regression_data,
    ):
        """
        Test generate_plot_mpl method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.liv_plot.LIVPlotter
        species_list : List of species to plot.
        nelements : int, Number of elements to include in plot.
        packets_mode : str, Packet mode, either 'virtual' or 'real'.
        xlog_scale : bool, If True, x-axis is scaled logarithmically.
        ylog_scale : bool, If True, y-axis is scaled logarithmically.
        num_bins : int, Number of bins for regrouping within the same range.
        velocity_range : tuple, Limits for the x-axis.
        """
        fig = plotter.generate_plot_mpl(
            species_list=species_list,
            nelements=nelements,
            packets_mode=packets_mode,
            xlog_scale=xlog_scale,
            ylog_scale=ylog_scale,
            num_bins=num_bins,
            velocity_range=velocity_range,
        )

        actual_species_name = plotter._species_name
        flat_list = [
            item for sublist in plotter._color_list for item in sublist
        ]
        actual_color_list = flat_list
        actual_step_x = list(plotter.step_x.value)
        actual_step_y = plotter.step_y

        actual = np.array(
            [str(item) for item in actual_species_name]
            + [str(item) for item in actual_color_list]
            + [str(item) for item in actual_step_x]
            + [str(item) for item in actual_step_y]
        )
        expected = regression_data.sync_ndarray(actual)
        npt.assert_array_equal(actual, expected)

    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"], None]
    )
    @pytest.mark.parametrize("nelements", [1, None])
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("xlog_scale", [True, False])
    @pytest.mark.parametrize("ylog_scale", [True, False])
    @pytest.mark.parametrize("num_bins", [10, 25])
    @pytest.mark.parametrize("velocity_range", [(18000, 25000)])
    def test_generate_plot_ply(
        self,
        plotter,
        species_list,
        nelements,
        packets_mode,
        xlog_scale,
        ylog_scale,
        num_bins,
        velocity_range,
        regression_data,
    ):
        """
        Test generate_plot_ply method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.liv_plot.LIVPlotter
        species_list : List of species to plot.
        nelements : int, Number of elements to include in plot.
        packets_mode : str, Packet mode, either 'virtual' or 'real'.
        xlog_scale : bool, If True, x-axis is scaled logarithmically.
        ylog_scale : bool, If True, y-axis is scaled logarithmically.
        num_bins : int, Number of bins for regrouping within the same range.
        velocity_range : tuple, Limits for the x-axis.
        """
        fig = plotter.generate_plot_ply(
            species_list=species_list,
            nelements=nelements,
            packets_mode=packets_mode,
            xlog_scale=xlog_scale,
            ylog_scale=ylog_scale,
            num_bins=num_bins,
            velocity_range=velocity_range,
        )

        actual_species_name = plotter._species_name
        flat_list = [
            item for sublist in plotter._color_list for item in sublist
        ]
        actual_color_list = flat_list
        actual_step_x = list(plotter.step_x.value)
        actual_step_y = plotter.step_y

        actual = np.array(
            [str(item) for item in actual_species_name]
            + [str(item) for item in actual_color_list]
            + [str(item) for item in actual_step_x]
            + [str(item) for item in actual_step_y]
        )
        expected = regression_data.sync_ndarray(actual)
        npt.assert_array_equal(actual, expected)
