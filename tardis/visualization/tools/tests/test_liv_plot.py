import numpy as np
import pytest
from numpy import testing as npt
from pandas import testing as pdt
from copy import deepcopy
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
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


def convert_to_native_type(obj):
    if isinstance(obj, dict):
        return {k: convert_to_native_type(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_native_type(i) for i in obj]
    elif isinstance(obj, np.int64):
        return int(obj)
    else:
        return obj


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
        request,
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
        subgroup_name = make_valid_name(request.node.callspec.id)
        plotter._parse_species_list(
            species_list=species_list,
            packets_mode=packets_mode,
            nelements=nelements,
        )
        regression_data_fname = (
            f"livplotter_parse_species_list_{subgroup_name}.h5"
        )

        expected = pd.read_hdf(regression_data_fname, "species_list")
        pdt.assert_frame_equal(plotter._species_list, expected)

        expected = pd.read_hdf(regression_data_fname, "keep_colour")
        pdt.assert_frame_equal(plotter._keep_colour, expected)

        expected = pd.read_hdf(regression_data_fname, "species_mapped")
        pdt.assert_frame_equal(
            convert_to_native_type(plotter._species_mapped), expected
        )

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"]]
    )
    @pytest.mark.parametrize("cmapname", ["jet"])
    @pytest.mark.parametrize("num_bins", [10, 25])
    @pytest.mark.parametrize("nelements", [1, None])
    def test_prepare_plot_data(
        self,
        request,
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
        subgroup_name = make_valid_name(request.node.callspec.id)
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
        plot_data_list = np.array(flat_list)
        regression_data_fname = (
            f"livplotter_prepare_plot_data_{subgroup_name}.h5"
        )

        expected = pd.read_hdf(regression_data_fname, "plot_data")
        pdt.assert_frame_equal(plot_data_list, expected)

        expected = pd.read_hdf(regression_data_fname, "plot_colors")
        pdt.assert_frame_equal(plotter.plot_colors, expected)

        expected = pd.read_hdf(regression_data_fname, "new_bin_edges")
        pdt.assert_frame_equal(plotter.new_bin_edges, expected)

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
        request,
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
        subgroup_name = make_valid_name("mpl" + request.node.callspec.id)
        fig = plotter.generate_plot_mpl(
            species_list=species_list,
            nelements=nelements,
            packets_mode=packets_mode,
            xlog_scale=xlog_scale,
            ylog_scale=ylog_scale,
            num_bins=num_bins,
            velocity_range=velocity_range,
        )
        fig_data = {
            "_species_name": plotter._species_name,
            "_color_list": plotter._color_list,
            "step_x": plotter.step_x,
            "step_y": plotter.step_y,
            "fig_data": [],
        }

        for index, data in enumerate(fig.get_children()):
            trace_data = {}
            if isinstance(data.get_label(), str):
                trace_data["label"] = data.get_label()
            if isinstance(data, PolyCollection):
                trace_data["paths"] = [
                    path.vertices for path in data.get_paths()
                ]
            if isinstance(data, Line2D):
                trace_data["xydata"] = data.get_xydata()
                trace_data["path"] = data.get_path().vertices
            fig_data["fig_data"].append(trace_data)

        regression_data_fname = (
            f"livplotter_generate_plot_mpl_{subgroup_name}.h5"
        )

        expected = pd.read_hdf(regression_data_fname, "fig_data")
        pdt.assert_frame_equal(fig_data, expected)

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
        request,
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
        subgroup_name = make_valid_name("ply" + request.node.callspec.id)
        fig = plotter.generate_plot_ply(
            species_list=species_list,
            nelements=nelements,
            packets_mode=packets_mode,
            xlog_scale=xlog_scale,
            ylog_scale=ylog_scale,
            num_bins=num_bins,
            velocity_range=velocity_range,
        )
        fig_data = {
            "_species_name": plotter._species_name,
            "_color_list": plotter._color_list,
            "step_x": plotter.step_x,
            "step_y": plotter.step_y,
            "fig_data": [],
        }

        for index, data in enumerate(fig.data):
            trace_data = {}
            if isinstance(data.name, str):
                trace_data["label"] = data.name
            if isinstance(data, go.Scatter):
                trace_data["x"] = data.x
                trace_data["y"] = data.y
            fig_data["fig_data"].append(trace_data)

        regression_data_fname = (
            f"livplotter_generate_plot_ply_{subgroup_name}.h5"
        )

        expected = pd.read_hdf(regression_data_fname, "fig_data")
        pdt.assert_frame_equal(fig_data, expected)
