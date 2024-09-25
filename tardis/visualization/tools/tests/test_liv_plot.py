from copy import deepcopy
from itertools import product

import astropy.units as u
import numpy as np
import pytest
from matplotlib.testing.compare import compare_images
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

from tardis.base import run_tardis
from tardis.io.util import HDFWriterMixin
from tardis.visualization.tools.liv_plot import LIVPlotter
from tardis.tests.fixtures.regression_data import RegressionData


class PlotDataHDF(HDFWriterMixin):
    """
    A class that writes plot data to HDF5 format using the HDFWriterMixin.
    """

    def __init__(self, **kwargs):
        """
        Initializes PlotDataHDF with arbitrary keyword arguments,
        storing them as attributes and adding their names to hdf_properties.

        Parameters:
        -----------
        **kwargs: Arbitrary keyword arguments representing properties to save.
        """
        self.hdf_properties = []
        for key, value in kwargs.items():
            setattr(self, key, value)
            self.hdf_properties.append(key)


@pytest.fixture(scope="module")
def simulation_simple(config_verysimple, atomic_dataset):
    """
    Fixture to create a simple TARDIS simulation.

    Parameters:
    -----------
    config_verysimple: A basic TARDIS configuration object.
    atomic_dataset: An atomic dataset to use in the simulation.

    Returns:
    --------
    A TARDIS simulation object.
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
    Fixture to create an LIVPlotter instance from a simulation.

    Parameters:
    -----------
    simulation_simple: A TARDIS simulation object.

    Returns:
    --------
    An LIVPlotter instance.
    """
    return LIVPlotter.from_simulation(simulation_simple)


class TestLIVPlotter:
    """Test the LIVPlotter class."""

    regression_data = None
    species_list = [["Si II", "Ca II", "C", "Fe I-V"], None]
    packet_wvl_range = [[3000, 9000] * u.AA]
    nelements = [1, None]
    packets_mode = ["virtual", "real"]
    num_bins = [10]
    velocity_range = [(18000, 25000)]
    cmapname = ["jet"]

    combinations = list(
        product(
            species_list,
            packet_wvl_range,
            packets_mode,
            nelements,
            num_bins,
            velocity_range,
            cmapname,
        )
    )

    @pytest.mark.parametrize(
        "attribute",
        [
            "_species_list",
            "_keep_colour",
            "_species_mapped",
        ],
    )
    def test_parse_species_list(
        self,
        request,
        plotter,
        attribute,
    ):
        """
        Test for the _parse_species_list method in LIVPlotter.

        Parameters:
        -----------
        request: Pytest's request fixture.
        plotter: The LIVPlotter instance.
        attribute: The attribute to test after parsing the species list.
        """
        regression_data = RegressionData(request)
        plotter._parse_species_list(
            packets_mode=self.packets_mode[0],
            species_list=self.species_list[0],
            nelements=self.nelements[0],
        )
        if attribute == "_species_mapped":
            plot_object = getattr(plotter, attribute)
            plot_object = [
                item
                for sublist in list(plot_object.values())
                for item in sublist
            ]
            data = regression_data.sync_ndarray(plot_object)
            np.testing.assert_allclose(plot_object, data)
        else:
            plot_object = getattr(plotter, attribute)
            data = regression_data.sync_ndarray(plot_object)
            np.testing.assert_allclose(plot_object, data)

    @pytest.fixture(scope="class", params=combinations)
    def plotter_prepare_plot_data(self, request, plotter):
        """
        Fixture to prepare plot data for a specific combination of parameters.

        Parameters:
        -----------
        request: Pytest's request fixture.
        plotter: The LIVPlotter instance.

        Returns:
        --------
        The plotter instance after preparing the plot data.
        """
        (
            species_list,
            packet_wvl_range,
            packets_mode,
            nelements,
            num_bins,
            _,
            cmapname,
        ) = request.param
        plotter._prepare_plot_data(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            species_list=species_list,
            cmapname=cmapname,
            num_bins=num_bins,
            nelements=nelements,
        )
        return plotter

    @pytest.mark.parametrize(
        "attribute",
        [
            "plot_data",
            "plot_colors",
            "new_bin_edges",
        ],
    )
    def test_prepare_plot_data(
        self,
        plotter_prepare_plot_data,
        request,
        attribute,
    ):
        """
        Test for the _prepare_plot_data method in LIVPlotter.

        Parameters:
        -----------
        plotter_prepare_plot_data: The plotter instance with prepared data.
        request: Pytest's request fixture.
        attribute: The attribute to test after preparing the plot data.
        """
        regression_data = RegressionData(request)
        if attribute == "plot_data" or attribute == "plot_colors":
            plot_object = getattr(plotter_prepare_plot_data, attribute)
            plot_object = [item for sublist in plot_object for item in sublist]
            if all(isinstance(item, u.Quantity) for item in plot_object):
                plot_object = [item.value for item in plot_object]
            data = regression_data.sync_ndarray(plot_object)
            np.testing.assert_allclose(plot_object, data)
        else:
            plot_object = getattr(plotter_prepare_plot_data, attribute)
            plot_object = plot_object.value
            data = regression_data.sync_ndarray(plot_object)
            np.testing.assert_allclose(plot_object, data)

    @pytest.fixture(scope="function", params=combinations)
    def plotter_generate_plot_mpl(self, request, plotter):
        """
        Fixture to generate a Matplotlib plot using the LIVPlotter.

        Parameters:
        -----------
        request: Pytest's request fixture.
        plotter: The LIVPlotter instance.

        Returns:
        --------
        A tuple containing the Matplotlib figure and the plotter instance.
        """
        (
            species_list,
            packet_wvl_range,
            packets_mode,
            nelements,
            num_bins,
            velocity_range,
            _,
        ) = request.param

        fig = plotter.generate_plot_mpl(
            species_list=species_list,
            packet_wvl_range=packet_wvl_range,
            nelements=nelements,
            packets_mode=packets_mode,
            num_bins=num_bins,
            velocity_range=velocity_range,
        )
        return fig, plotter

    @pytest.fixture(scope="function")
    def generate_plot_mpl_hdf(self, plotter_generate_plot_mpl):
        """
        Fixture to generate and store plot data for Matplotlib in HDF5 format.

        Parameters:
        -----------
        plotter_generate_plot_mpl: The Matplotlib plotter fixture.

        Returns:
        --------
        A PlotDataHDF instance containing the plot data.
        """
        fig, plotter = plotter_generate_plot_mpl

        color_list = [
            item for subitem in plotter._color_list for item in subitem
        ]
        property_group = {
            "_species_name": plotter._species_name,
            "_color_list": color_list,
            "step_x": plotter.step_x.value,
            "step_y": plotter.step_y,
        }
        for index1, data in enumerate(fig.get_children()):
            if isinstance(data.get_label(), str):
                property_group[
                    "label" + str(index1)
                ] = data.get_label().encode()
            if isinstance(data, Line2D):
                property_group["data" + str(index1)] = data.get_xydata()
                property_group[
                    "linepath" + str(index1)
                ] = data.get_path().vertices
            if isinstance(data, PolyCollection):
                for index2, path in enumerate(data.get_paths()):
                    property_group[
                        "polypath" + "ind_" + str(index1) + "ind_" + str(index2)
                    ] = path.vertices

        plot_data = PlotDataHDF(**property_group)
        return plot_data

    def test_generate_plot_mpl(
        self, generate_plot_mpl_hdf, plotter_generate_plot_mpl, request
    ):
        """
        Test for the generate_plot_mpl method in LIVPlotter.

        Compares generated plot data with regression data.

        Parameters:
        -----------
        generate_plot_mpl_hdf: The PlotDataHDF fixture for Matplotlib.
        plotter_generate_plot_mpl: The Matplotlib plotter fixture.
        request: Pytest's request fixture.
        """
        fig, _ = plotter_generate_plot_mpl
        regression_data = RegressionData(request)
        expected = regression_data.sync_hdf_store(generate_plot_mpl_hdf)
        for item in ["_species_name", "_color_list", "step_x", "step_y"]:
            expected_values = expected.get(
                "plot_data_hdf/" + item
            ).values.flatten()
            actual_values = getattr(generate_plot_mpl_hdf, item)

            if np.issubdtype(expected_values.dtype, np.number):
                np.testing.assert_allclose(
                    expected_values,
                    actual_values,
                    rtol=0.3,
                    atol=3,
                )
            else:
                assert np.array_equal(expected_values, actual_values)

        labels = expected["plot_data_hdf/scalars"]
        for index1, data in enumerate(fig.get_children()):
            if isinstance(data.get_label(), str):
                assert (
                    getattr(labels, "label" + str(index1)).decode()
                    == data.get_label()
                )
            if isinstance(data, Line2D):
                np.testing.assert_allclose(
                    data.get_xydata(),
                    expected.get("plot_data_hdf/" + "data" + str(index1)),
                    rtol=0.3,
                    atol=3,
                )
                np.testing.assert_allclose(
                    data.get_path().vertices,
                    expected.get("plot_data_hdf/" + "linepath" + str(index1)),
                    rtol=1,
                    atol=3,
                )
            if isinstance(data, PolyCollection):
                for index2, path in enumerate(data.get_paths()):
                    np.testing.assert_almost_equal(
                        path.vertices,
                        expected.get(
                            "plot_data_hdf/"
                            + "polypath"
                            + "ind_"
                            + str(index1)
                            + "ind_"
                            + str(index2)
                        ),
                    )

    def test_mpl_image(self, plotter_generate_plot_mpl, tmp_path, request):
        """
        Test to compare the generated Matplotlib images with the expected ones.

        Parameters:
        -----------
        plotter_generate_plot_mpl: The Matplotlib plotter fixture.
        request: Pytest's request fixture.
        recwarn: Pytest's warning recording fixture.
        pytestconfig: Pytest's configuration fixture.
        """
        regression_data = RegressionData(request)
        fig, _ = plotter_generate_plot_mpl
        regression_data.fpath.parent.mkdir(parents=True, exist_ok=True)
        fig.figure.savefig(tmp_path / f"{regression_data.fname_prefix}.png")

        if regression_data.enable_generate_reference:
            fig.figure.savefig(
                regression_data.absolute_regression_data_dir
                / f"{regression_data.fname_prefix}.png"
            )
            pytest.skip("Skipping test to generate reference data")
        else:
            expected = str(
                regression_data.absolute_regression_data_dir
                / f"{regression_data.fname_prefix}.png"
            )
            actual = str(tmp_path / f"{regression_data.fname_prefix}.png")
            compare_images(expected, actual, tol=0.001)

    @pytest.fixture(scope="function", params=combinations)
    def plotter_generate_plot_ply(self, request, plotter):
        """
        Fixture to generate a Plotly plot using the LIVPlotter.

        Parameters:
        -----------
        request: Pytest's request fixture.
        plotter: The LIVPlotter instance.

        Returns:
        --------
        A tuple containing the Plotly figure and the plotter instance.
        """
        (
            species_list,
            packet_wvl_range,
            packets_mode,
            nelements,
            num_bins,
            velocity_range,
            _,
        ) = request.param

        fig = plotter.generate_plot_ply(
            species_list=species_list,
            packet_wvl_range=packet_wvl_range,
            nelements=nelements,
            packets_mode=packets_mode,
            num_bins=num_bins,
            velocity_range=velocity_range,
        )
        return fig, plotter

    @pytest.fixture(scope="function")
    def generate_plot_plotly_hdf(self, plotter_generate_plot_ply):
        """
        Fixture to generate and store plot data for Matplotlib in HDF5 format.

        Parameters:
        -----------
        plotter_generate_plot_ply: The Plotly plotter fixture.

        Returns:
        --------
        A PlotDataHDF instance containing the plot data.
        """
        fig, plotter = plotter_generate_plot_ply

        color_list = [
            item for subitem in plotter._color_list for item in subitem
        ]
        property_group = {
            "_species_name": plotter._species_name,
            "_color_list": color_list,
            "step_x": plotter.step_x.value,
            "step_y": plotter.step_y,
        }
        for index, data in enumerate(fig.data):
            group = "_" + str(index)
            if data.stackgroup:
                property_group[group + "stackgroup"] = data.stackgroup.encode()
            if data.name:
                property_group[group + "name"] = data.name.encode()
            property_group[group + "x"] = data.x
            property_group[group + "y"] = data.y
        plot_data = PlotDataHDF(**property_group)
        return plot_data

    def test_generate_plot_ply(
        self, generate_plot_plotly_hdf, plotter_generate_plot_ply, request
    ):
        """
        Test for the generate_plot_mpl method in LIVPlotter.

        Compares generated plot data with regression data.

        Parameters:
        ----------
        generate_plot_plotly_hdf: The PlotDataHDF fixture for Plotly.
        plotter_generate_plot_mpl: The Plotly plotter fixture.
        request: Pytest's request fixture.
        """
        fig, _ = plotter_generate_plot_ply
        regression_data = RegressionData(request)
        expected = regression_data.sync_hdf_store(generate_plot_plotly_hdf)

        for item in ["_species_name", "_color_list", "step_x", "step_y"]:
            expected_values = expected.get(
                "plot_data_hdf/" + item
            ).values.flatten()
            actual_values = getattr(generate_plot_plotly_hdf, item)

            if np.issubdtype(expected_values.dtype, np.number):
                np.testing.assert_allclose(
                    expected_values,
                    actual_values,
                    rtol=0.3,
                    atol=3,
                )
            else:
                assert np.array_equal(expected_values, actual_values)
        for index, data in enumerate(fig.data):
            group = "plot_data_hdf/" + "_" + str(index)
            if data.stackgroup:
                assert (
                    data.stackgroup
                    == getattr(
                        expected["/plot_data_hdf/scalars"],
                        "_" + str(index) + "stackgroup",
                    ).decode()
                )
            if data.name:
                assert (
                    data.name
                    == getattr(
                        expected["/plot_data_hdf/scalars"],
                        "_" + str(index) + "name",
                    ).decode()
                )
            np.testing.assert_allclose(
                data.x,
                expected.get(group + "x").values.flatten(),
                rtol=0.3,
                atol=3,
            )
            np.testing.assert_allclose(
                data.y,
                expected.get(group + "y").values.flatten(),
                rtol=0.3,
                atol=3,
            )
