from itertools import product
from copy import deepcopy

import astropy.units as u
import numpy as np
import pytest
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.testing.compare import compare_images

import pandas as pd

from tardis.visualization.tools.liv_plot import LIVPlotter
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow
from tardisbase.testing.regression_data.regression_data import PlotDataHDF

RELATIVE_TOLERANCE_LIV=1e-12

@pytest.fixture(scope="class")
def plotter(simulation_simple_tracked):
    """
    Fixture to create an LIVPlotter instance from a simulation.

    Parameters:
    -----------
    simulation_simple_tracked: A TARDIS simulation object.

    Returns:
    --------
    An LIVPlotter instance.
    """
    return LIVPlotter.from_simulation(simulation_simple_tracked)


class TestLIVPlotter:
    """Test the LIVPlotter class."""

    regression_data = None
    species_list = [["Si II", "Ca II", "C", "Fe I-V"], None]
    packet_wvl_range = [[3000, 9000] * u.AA]
    nelements = [1, None]
    num_bins = [10]
    velocity_range = [(18000, 25000)]
    cmapname = ["jet"]

    combinations = list(
        product(
            species_list,
            packet_wvl_range,
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
        regression_data,
        plotter,
        attribute,
    ):
        """
        Test for the _parse_species_list method in LIVPlotter.

        Parameters:
        -----------
        regression_data : _pytest.fixtures.RegressionData
        plotter: The LIVPlotter instance.
        attribute: The attribute to test after parsing the species list.
        """
        plotter._parse_species_list(
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
            np.testing.assert_allclose(plot_object, data, atol=0, rtol=RELATIVE_TOLERANCE_LIV)
        else:
            plot_object = getattr(plotter, attribute)
            data = regression_data.sync_ndarray(plot_object)
            np.testing.assert_allclose(plot_object, data, atol=0, rtol=RELATIVE_TOLERANCE_LIV)

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
            nelements,
            num_bins,
            _,
            cmapname,
        ) = request.param
        plotter._prepare_plot_data(
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
        regression_data,
        attribute,
    ):
        """
        Test for the _prepare_plot_data method in LIVPlotter.

        Parameters:
        -----------
        plotter_prepare_plot_data: The plotter instance with prepared data.
        regression_data: _pytest.fixtures.RegressionData
        attribute: The attribute to test after preparing the plot data.
        """
        if attribute == "plot_data" or attribute == "plot_colors":
            plot_object = getattr(plotter_prepare_plot_data, attribute)
            plot_object = [item for sublist in plot_object for item in sublist]
            if all(isinstance(item, u.Quantity) for item in plot_object):
                plot_object = [item.value for item in plot_object]
            data = regression_data.sync_ndarray(plot_object)
            np.testing.assert_allclose(plot_object, data, atol=0, rtol=RELATIVE_TOLERANCE_LIV)
        else:
            plot_object = getattr(plotter_prepare_plot_data, attribute)
            plot_object = plot_object.value
            data = regression_data.sync_ndarray(plot_object)
            np.testing.assert_allclose(plot_object, data, atol=0, rtol=RELATIVE_TOLERANCE_LIV)

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
            nelements,
            num_bins,
            velocity_range,
            _,
        ) = request.param

        fig = plotter.generate_plot_mpl(
            species_list=species_list,
            packet_wvl_range=packet_wvl_range,
            nelements=nelements,
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
        self, generate_plot_mpl_hdf, plotter_generate_plot_mpl, regression_data
    ):
        """
        Test for the generate_plot_mpl method in LIVPlotter.

        Compares generated plot data with regression data.

        Parameters:
        -----------
        generate_plot_mpl_hdf: The PlotDataHDF fixture for Matplotlib.
        plotter_generate_plot_mpl: The Matplotlib plotter fixture.
        regression_data : _pytest.fixtures.RegressionData
        """
        fig, _ = plotter_generate_plot_mpl
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
                    atol=0, rtol=RELATIVE_TOLERANCE_LIV
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
                    atol=0, rtol=RELATIVE_TOLERANCE_LIV
                )
                np.testing.assert_allclose(
                    data.get_path().vertices,
                    expected.get("plot_data_hdf/" + "linepath" + str(index1)),
                    atol=0, rtol=RELATIVE_TOLERANCE_LIV
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
        expected.close()

    def test_mpl_image(self, plotter_generate_plot_mpl, tmp_path, regression_data):
        """
        Test to compare the generated Matplotlib images with the expected ones.

        Parameters:
        -----------
        plotter_generate_plot_mpl : Fixture that returns a Matplotlib figure and axes (fig, ax).
        tmp_path : Temporary directory provided by pytest for saving generated output
        regression_data : _pytest.fixtures.RegressionData.
        """
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
            compare_images(expected, actual, tol=1e-3)

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
            nelements,
            num_bins,
            velocity_range,
            _,
        ) = request.param

        fig = plotter.generate_plot_ply(
            species_list=species_list,
            packet_wvl_range=packet_wvl_range,
            nelements=nelements,
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
        self, generate_plot_plotly_hdf, plotter_generate_plot_ply, regression_data
    ):
        """
        Test for the generate_plot_mpl method in LIVPlotter.

        Compares generated plot data with regression data.

        Parameters:
        ----------
        generate_plot_plotly_hdf: The PlotDataHDF fixture for Plotly.
        plotter_generate_plot_mpl: The Plotly plotter fixture.
        regression_data : _pytest.fixtures.RegressionData.
        """
        fig, _ = plotter_generate_plot_ply
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
                    atol=0, rtol=RELATIVE_TOLERANCE_LIV
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
                atol=0, rtol=RELATIVE_TOLERANCE_LIV
            )
            np.testing.assert_allclose(
                data.y,
                expected.get(group + "y").values.flatten(),
                atol=0, rtol=RELATIVE_TOLERANCE_LIV
            )
        expected.close()

    @pytest.fixture(scope="class")
    def plotter_from_workflow(self, workflow_simple):
        return LIVPlotter.from_workflow(workflow_simple)

    @pytest.fixture(scope="class")
    def liv_regression_data(self, tardis_regression_path):
        return tardis_regression_path / "tardis/visualization/tools/tests/test_liv_plot/test_liv_plotter"

    @pytest.fixture(scope="class", params=list(enumerate(combinations)))
    def plotter_prepare_plot_data_from_workflow(self, request, plotter_from_workflow):
        param_idx, param = request.param
        (
            species_list,
            packet_wvl_range,
            nelements,
            num_bins,
            _,
            cmapname,
        ) = param
        plotter_from_workflow._prepare_plot_data(
            packet_wvl_range=packet_wvl_range,
            species_list=species_list,
            cmapname=cmapname,
            num_bins=num_bins,
            nelements=nelements,
        )
        plotter_from_workflow._param_idx = param_idx
        return plotter_from_workflow

    def test_prepare_plot_data_workflow_vs_regression(
        self, plotter_prepare_plot_data_from_workflow, liv_regression_data
    ):
        param_idx = plotter_prepare_plot_data_from_workflow._param_idx
        
        for attribute in ["plot_data", "plot_colors", "new_bin_edges"]:
            if attribute == "plot_data" or attribute == "plot_colors":
                regression_file = liv_regression_data / f"test_prepare_plot_data__plotter_prepare_plot_data{param_idx}-{attribute}__.npy"
                expected = np.load(regression_file)
                
                plot_object = getattr(plotter_prepare_plot_data_from_workflow, attribute)
                plot_object = [item for sublist in plot_object for item in sublist]
                if all(isinstance(item, u.Quantity) for item in plot_object):
                    plot_object = [item.value for item in plot_object]
                    
                np.testing.assert_allclose(plot_object, expected, atol=0, rtol=RELATIVE_TOLERANCE_LIV)
            else:
                regression_file = liv_regression_data / f"test_prepare_plot_data__plotter_prepare_plot_data{param_idx}-{attribute}__.npy"
                expected = np.load(regression_file)
                
                plot_object = getattr(plotter_prepare_plot_data_from_workflow, attribute)
                plot_object = plot_object.value
                
                np.testing.assert_allclose(plot_object, expected, atol=0, rtol=RELATIVE_TOLERANCE_LIV)

    @pytest.fixture(scope="class", params=list(enumerate(combinations)))
    def plotter_generate_plot_mpl_from_workflow(self, request, plotter_from_workflow):
        param_idx, param = request.param
        (
            species_list,
            packet_wvl_range,
            nelements,
            num_bins,
            velocity_range,
            _,
        ) = param

        fig = plotter_from_workflow.generate_plot_mpl(
            species_list=species_list,
            packet_wvl_range=packet_wvl_range,
            nelements=nelements,
            num_bins=num_bins,
            velocity_range=velocity_range,
        )
        plotter_from_workflow._param_idx = param_idx
        return fig, plotter_from_workflow

    def test_generate_plot_mpl_workflow_vs_regression(
        self, plotter_generate_plot_mpl_from_workflow, liv_regression_data
    ):
        _, plotter = plotter_generate_plot_mpl_from_workflow
        param_idx = plotter._param_idx
        regression_file = liv_regression_data / f"test_generate_plot_mpl__plotter_generate_plot_mpl{param_idx}__.h5"
        
        # Compare species names and color lists
        expected_species = pd.read_hdf(regression_file, key="plot_data_hdf/_species_name")
        expected_colors = pd.read_hdf(regression_file, key="plot_data_hdf/_color_list")
        expected_step_x = pd.read_hdf(regression_file, key="plot_data_hdf/step_x")
        expected_step_y = pd.read_hdf(regression_file, key="plot_data_hdf/step_y")
        
        np.testing.assert_array_equal(plotter._species_name, expected_species.values.flatten())
        
        color_list = [item for subitem in plotter._color_list for item in subitem]
        np.testing.assert_allclose(color_list, expected_colors.values.flatten(), atol=0, rtol=RELATIVE_TOLERANCE_LIV)
        
        np.testing.assert_allclose(plotter.step_x.value, expected_step_x.values.flatten(), atol=0, rtol=RELATIVE_TOLERANCE_LIV)
        np.testing.assert_allclose(plotter.step_y, expected_step_y.values.flatten(), atol=0, rtol=RELATIVE_TOLERANCE_LIV)

    def test_workflow_simulation_data_identical(self, plotter, plotter_from_workflow):
        # Calculate plotting data with identical parameters
        params = {
            "packet_wvl_range": [3000, 9000] * u.AA,
            "species_list": None,
            "cmapname": "jet",
            "num_bins": None,
            "nelements": None
        }

        plotter._prepare_plot_data(**params)
        plotter_from_workflow._prepare_plot_data(**params)
        
        # Test that plot data arrays are identical
        for i, (sim_data, workflow_data) in enumerate(zip(plotter.plot_data, plotter_from_workflow.plot_data)):
            np.testing.assert_allclose(
                [item.value for item in sim_data],
                [item.value for item in workflow_data],
                atol=0, rtol=RELATIVE_TOLERANCE_LIV
            )

    def test_mpl_image_workflow(self, plotter_generate_plot_mpl_from_workflow, tmp_path, liv_regression_data):
        fig, plotter = plotter_generate_plot_mpl_from_workflow
        param_idx = plotter._param_idx
        
        # Save actual image
        actual_image_path = tmp_path / f"test_mpl_image_workflow_{param_idx}.png"
        fig.figure.savefig(actual_image_path)
        
        # Path to expected image
        expected_image_path = liv_regression_data / f"test_mpl_image__plotter_generate_plot_mpl{param_idx}__.png"
        
        # Compare images
        compare_images(str(expected_image_path), str(actual_image_path), tol=1e-3)
