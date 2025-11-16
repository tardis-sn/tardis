"""Tests for SDEC Plots."""
from itertools import product
from copy import deepcopy

import astropy
import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.testing.compare import compare_images

from tardisbase.testing.regression_data.regression_data import PlotDataHDF
from tardis.visualization.tools.sdec_plot import SDECPlotter
from tardis.workflows.standard_tardis_workflow import StandardTARDISWorkflow

RELATIVE_TOLERANCE_SDEC=1e-12

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
    testid = testid.replace("-", "_")
    testid = "_" + testid
    return testid


@pytest.fixture()
def sdec_regression_data(tardis_regression_path):
    # workflow tests for the SDEC plot use the existing regression data from the 
    # Simulation since both objects produce the same plot with same config.
    return tardis_regression_path / "tardis/visualization/tools/tests/test_sdec_plot/test_sdec_plotter"


class TestSDECPlotter:
    """Test the SDECPlotter class."""

    regression_data = None
    distance = [10 * u.Mpc, None]
    packet_wvl_range = [[500, 9000] * u.AA]
    species_list = [["Si II", "Ca II", "C", "Fe I-V"]]
    packets_mode = ["real", "virtual"]
    nelements = [1, None]
    show_modeled_spectrum = [True, False]

    combinations = list(
        product(
            distance,
            packet_wvl_range,
            species_list,
            packets_mode,
            nelements,
            show_modeled_spectrum,
        )
    )

    plotting_data_attributes = {
        "attributes_np": [
            "plot_frequency_bins",
            "plot_wavelength",
            "plot_frequency",
            "modeled_spectrum_luminosity",
            "packet_wvl_range_mask",
            "emission_species",
            "absorption_species",
        ],
        "attributes_df": [
            "absorption_luminosities_df",
            "emission_luminosities_df",
            "total_luminosities_df",
        ],
    }
    plotting_data_attributes = [
        [key, value]
        for key, values in plotting_data_attributes.items()
        for value in values
    ]

    @pytest.fixture(scope="class")
    def plotter(self, simulation_simple_tracked):
        """
        Create a SDECPlotter object.

        Parameters
        ----------
        simulation_simple_tracked : tardis.simulation.base.Simulation
            Simulation object.

        Returns
        -------
        tardis.visualization.tools.sdec_plot.SDECPlotter
        """
        return SDECPlotter.from_simulation(simulation_simple_tracked)

    @pytest.fixture(scope="class")
    def observed_spectrum(self):
        """
        Return the observed spectrum.

        Returns
        -------
        Tuple of two astropy.units.quantity.Quantity values.
        """
        test_data = np.loadtxt(
            "tardis/visualization/tools/tests/data/observed_spectrum_test_data.dat"
        )
        observed_spectrum_wavelength, observed_spectrum_flux = test_data.T
        observed_spectrum_wavelength = observed_spectrum_wavelength * u.AA
        observed_spectrum_flux = (
            observed_spectrum_flux * u.erg / (u.s * u.cm**2 * u.AA)
        )
        return observed_spectrum_wavelength, observed_spectrum_flux

    @pytest.mark.parametrize(
        "attribute", ["_full_species_list", "_species_list", "_keep_colour"]
    )
    def test_parse_species_list(self, regression_data, plotter, attribute):
        """
        Test _parse_species_list method.

        Parameters
        ----------
        regression_data : _pytest.fixtures.RegressionData
        plotter : tardis.visualization.tools.sdec_plot.SDECPlotter
        attribute : parameter to be tested
        """
        # THIS NEEDS TO BE RUN FIRST. NOT INDEPENDENT TESTS
        plotter._parse_species_list(self.species_list[0])
        data = regression_data.sync_ndarray(getattr(plotter, attribute))
        if attribute == "_full_species_list":
            np.testing.assert_equal(getattr(plotter, attribute), data)
        else:
            np.testing.assert_allclose(getattr(plotter, attribute), data, atol=0, rtol=RELATIVE_TOLERANCE_SDEC)

    @pytest.fixture(scope="class", params=combinations)
    def plotter_calculate_plotting_data(self, request, plotter):
        (
            distance,
            packet_wvl_range,
            _,
            packets_mode,
            nelements,
            _,
        ) = request.param
        # plotter._parse_species_list(species_list=None)
        if packets_mode == "virtual":
          pytest.skip("Skipping tests for virtual packets mode")

        # the tests connected to this fixture wont run individually, 
        # since this needs the parse species list to run

        plotter._calculate_plotting_data(
            packets_mode, packet_wvl_range, distance, nelements
        )
        return plotter

    @pytest.fixture(scope="class")
    def calculate_plotting_data_hdf(
        self, plotter_calculate_plotting_data
    ):
        property_group = {}
        for _, attribute_name in self.plotting_data_attributes:
            plot_object = getattr(
                plotter_calculate_plotting_data, attribute_name
            )
            property_group[attribute_name] = plot_object
        plot_data = PlotDataHDF(**property_group)
        return plot_data

    def test_calculate_plotting_data(
        self,
        plotter_calculate_plotting_data,
        calculate_plotting_data_hdf,
        regression_data,
    ):
        expected = regression_data.sync_hdf_store(calculate_plotting_data_hdf)
        group = "plot_data_hdf/"
        for attribute_type, attribute_name in self.plotting_data_attributes:
            plot_object = getattr(
                plotter_calculate_plotting_data, attribute_name
            )
            if attribute_type == "attributes_np":
                if isinstance(plot_object, astropy.units.quantity.Quantity):
                    plot_object = plot_object.cgs.value
                np.testing.assert_allclose(
                    plot_object, expected.get(group + attribute_name), atol=0, rtol=RELATIVE_TOLERANCE_SDEC
                )
            if attribute_type == "attributes_pd":
                pd.testing.assert_frame_equal(
                    plot_object, expected.get(group + attribute_name), atol=0, rtol=RELATIVE_TOLERANCE_SDEC
                )
        expected.close()

    @pytest.fixture(scope="class", params=list(enumerate(combinations)))
    def plotter_generate_plot_mpl(self, request, observed_spectrum, plotter, simulation_simple_tracked):
        param_idx, param = request.param
        
        (
            distance,
            packet_wvl_range,
            species_list,
            packets_mode,
            nelements,
            show_modeled_spectrum,
        ) = param

        if distance is None:
            observed_spectrum = None
        
        # plotter = SDECPlotter.from_simulation(simulation_simple_tracked)
        if packets_mode == "virtual":
          pytest.skip("Skipping tests for virtual packets mode")


        fig = plotter.generate_plot_mpl(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            show_modeled_spectrum=show_modeled_spectrum,
            observed_spectrum=observed_spectrum,
            nelements=nelements,
            species_list=species_list,
        )
        plotter._param_idx = param_idx
        return fig, plotter

    @pytest.fixture(scope="class")
    def generate_plot_mpl_hdf(self, plotter_generate_plot_mpl):
        fig, plotter = plotter_generate_plot_mpl

        color_list = [
            item for subitem in plotter._color_list for item in subitem
        ]
        property_group = {
            "_species_name": plotter._species_name,
            "_color_list": color_list,
        }
        for index1, data in enumerate(fig.get_children()):
            if isinstance(data.get_label(), str):
                property_group[
                    "label" + str(index1)
                ] = data.get_label().encode()
            # save line plots
            if isinstance(data, Line2D):
                property_group["data" + str(index1)] = data.get_xydata()
                property_group[
                    "linepath" + str(index1)
                ] = data.get_path().vertices

            # save artists which correspond to element contributions
            if isinstance(data, PolyCollection):
                for index2, path in enumerate(data.get_paths()):
                    property_group[
                        "polypath" + "ind_" + str(index1) + "ind_" + str(index2)
                    ] = path.vertices

        plot_data = PlotDataHDF(**property_group)
        return plot_data

    def test_generate_plot_mpl(
        self, generate_plot_mpl_hdf, plotter_generate_plot_mpl, regression_data, sdec_regression_data
    ):
        fig, plotter = plotter_generate_plot_mpl
        param_idx = plotter._param_idx
        regression_file = f"test_generate_plot_mpl__plotter_generate_plot_ply{param_idx}__.h5"
        
        regression_data.fname = regression_file
        # expected = regression_data.sync_hdf_store(generate_plot_plotly_hdf)
        expected = pd.HDFStore(sdec_regression_data / regression_file, mode='r')

        # expected = regression_data.sync_hdf_store(generate_plot_mpl_hdf)

        for item in ["_species_name", "_color_list"]:
            np.testing.assert_array_equal(
                expected.get("plot_data_hdf/" + item).values.flatten(),
                getattr(generate_plot_mpl_hdf, item),
            )

        labels = expected["plot_data_hdf/scalars"]
        for index1, data in enumerate(fig.get_children()):
            if isinstance(data.get_label(), str):
                assert (
                    getattr(labels, "label" + str(index1)).decode()
                    == data.get_label()
                )
            # save line plots
            if isinstance(data, Line2D):
                np.testing.assert_allclose(
                    data.get_xydata(),
                    expected.get("plot_data_hdf/" + "data" + str(index1)),
                    atol=0, rtol=RELATIVE_TOLERANCE_SDEC
                )
                np.testing.assert_allclose(
                    data.get_path().vertices,
                    expected.get("plot_data_hdf/" + "linepath" + str(index1)),
                    atol=0, rtol=RELATIVE_TOLERANCE_SDEC
                )
            # save artists which correspond to element contributions
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

    @pytest.fixture(scope="class", params=list(enumerate(combinations)))
    def plotter_generate_plot_ply(self, request, observed_spectrum, plotter):
        param_idx, param = request.param

        (
            distance,
            packet_wvl_range,
            species_list,
            packets_mode,
            nelements,
            show_modeled_spectrum,
        ) = param
        if distance is None:
            observed_spectrum = None

        if packets_mode == "virtual":
          pytest.skip("Skipping tests for virtual packets mode")

        fig = plotter.generate_plot_ply(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            show_modeled_spectrum=show_modeled_spectrum,
            observed_spectrum=observed_spectrum,
            nelements=nelements,
            species_list=species_list,
        )
        # import pdb; pdb.set_trace()
        plotter._param_idx = param_idx
        return fig, plotter

    @pytest.fixture(scope="class")
    def generate_plot_plotly_hdf(self, plotter_generate_plot_ply):
        fig, plotter = plotter_generate_plot_ply

        color_list = [
            item for subitem in plotter._color_list for item in subitem
        ]
        property_group = {
            "_species_name": plotter._species_name,
            "_color_list": color_list,
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
        self, generate_plot_plotly_hdf, plotter_generate_plot_ply, regression_data, sdec_regression_data
    ):
        fig, plotter = plotter_generate_plot_ply
        param_idx = plotter._param_idx
        regression_file = f"test_generate_plot_mpl__plotter_generate_plot_ply{param_idx}__.h5"
        
        regression_data.fname = regression_file
        # expected = regression_data.sync_hdf_store(generate_plot_plotly_hdf)
        expected = pd.HDFStore(sdec_regression_data / regression_file, mode='r')

        for item in ["_species_name", "_color_list"]:
            np.testing.assert_array_equal(
                expected.get("plot_data_hdf/" + item).values.flatten(),
                getattr(generate_plot_plotly_hdf, item),
            )

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
                data.x, expected.get(group + "x").values.flatten(), atol=0, rtol=RELATIVE_TOLERANCE_SDEC
            )
            np.testing.assert_allclose(
                data.y, expected.get(group + "y").values.flatten(), atol=0, rtol=RELATIVE_TOLERANCE_SDEC
            )

        expected.close()

    def test_mpl_image(self, plotter_generate_plot_mpl, tmp_path, regression_data):
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

    def test_make_colorbar_labels(self, plotter):
        plotter._parse_species_list(None)
        plotter._calculate_plotting_data(
            packets_mode="virtual",
            packet_wvl_range=[500, 9000] * u.AA,
            distance=None,
            nelements=None,
        )
        plotter._make_colorbar_labels()
        assert isinstance(plotter._species_name, list)
        assert all(isinstance(label, str) for label in plotter._species_name)

    @pytest.fixture(scope="class")
    def plotter_from_workflow(self, workflow_simple_tracked):
        return SDECPlotter.from_workflow(workflow_simple_tracked)

    def test_from_workflow_vs_from_simulation_data_consistency(
        self, plotter, plotter_from_workflow
    ):
        # Test key attributes are equal
        np.testing.assert_allclose(
            plotter.t_inner.value, plotter_from_workflow.t_inner.value, atol=0, rtol=RELATIVE_TOLERANCE_SDEC
        )
        np.testing.assert_allclose(
            plotter.r_inner.value, plotter_from_workflow.r_inner.value, atol=0, rtol=RELATIVE_TOLERANCE_SDEC
        )
        np.testing.assert_allclose(
            plotter.time_of_simulation.value, 
            plotter_from_workflow.time_of_simulation.value, 
            atol=0, rtol=RELATIVE_TOLERANCE_SDEC
        )
        
        # Test packet data structures exist for both modes
        for mode in ["real", "virtual"]:
            assert plotter.packet_data[mode]["packets_df"] is not None
            assert plotter_from_workflow.packet_data[mode]["packets_df"] is not None
            assert plotter.spectrum[mode] is not None
            assert plotter_from_workflow.spectrum[mode] is not None


    def test_from_workflow_method_functionality(self, plotter_from_workflow):
        plotter_from_workflow._parse_species_list(None)
        plotter_from_workflow._calculate_plotting_data(
            packets_mode="virtual",
            packet_wvl_range=[500, 9000] * u.AA, 
            distance=None,
            nelements=None
        )
        
        # Test that basic plotting data exists and has reasonable values
        assert plotter_from_workflow.emission_luminosities_df is not None
        assert plotter_from_workflow.absorption_luminosities_df is not None
        assert len(plotter_from_workflow.emission_luminosities_df) > 0
        assert len(plotter_from_workflow.absorption_luminosities_df) > 0
        
        # Test that matplotlib plot can be generated
        fig = plotter_from_workflow.generate_plot_mpl(packets_mode="virtual")
        assert fig is not None

    @pytest.fixture(scope="class", params=list(enumerate(combinations)))
    def plotter_calculate_plotting_data_from_workflow(self, request, plotter_from_workflow):
        param_idx, param = request.param
        (
            distance,
            packet_wvl_range,
            species_list,
            packets_mode,
            nelements,
            _,
        ) = param
        if packets_mode == "virtual":
          pytest.skip("Skipping tests for virtual packets mode")

        # we need to parse this
        plotter_from_workflow._parse_species_list(species_list)
        
        plotter_from_workflow._calculate_plotting_data(
            packets_mode, packet_wvl_range, distance, nelements
        )
        plotter_from_workflow._param_idx = param_idx
        return plotter_from_workflow


    def test_calculate_plotting_data_workflow_vs_regression(
        self, plotter_calculate_plotting_data_from_workflow, sdec_regression_data
    ):
        param_idx = plotter_calculate_plotting_data_from_workflow._param_idx
        regression_file = sdec_regression_data / f"test_calculate_plotting_data__plotter_calculate_plotting_data{param_idx}__.h5"
        
        for attribute_type, attribute_name in self.plotting_data_attributes:
            plot_object = getattr(plotter_calculate_plotting_data_from_workflow, attribute_name)
            if attribute_type == "attributes_np":
                expected = pd.read_hdf(regression_file, key=f"plot_data_hdf/{attribute_name}", mode='r')
                if isinstance(plot_object, astropy.units.quantity.Quantity):
                    plot_object = plot_object.cgs.value
                # Handle array shape differences
                if plot_object.ndim > 1:
                    plot_object = plot_object.flatten()
                np.testing.assert_allclose(plot_object, expected.values.flatten(), atol=0, rtol=RELATIVE_TOLERANCE_SDEC)
            elif attribute_type == "attributes_df":
                expected_df = pd.read_hdf(regression_file, key=f"plot_data_hdf/{attribute_name}", mode='r')
                pd.testing.assert_frame_equal(plot_object, expected_df, atol=0, rtol=RELATIVE_TOLERANCE_SDEC)

    @pytest.fixture(scope="class", params=list(enumerate(combinations)))
    def plotter_generate_plot_mpl_from_workflow(self, request, observed_spectrum, plotter_from_workflow):
        param_idx, param = request.param
        (
            distance,
            packet_wvl_range,
            species_list,
            packets_mode,
            nelements,
            show_modeled_spectrum,
        ) = param
        if distance is None:
            observed_spectrum = None
        
        if packets_mode == "virtual":
          pytest.skip("Skipping tests for virtual packets mode")


        fig = plotter_from_workflow.generate_plot_mpl(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            show_modeled_spectrum=show_modeled_spectrum,
            observed_spectrum=observed_spectrum,
            nelements=nelements,
            species_list=species_list,
        )
        plotter_from_workflow._param_idx = param_idx
        return fig, plotter_from_workflow

    def test_generate_plot_mpl_workflow_vs_regression(
        self, plotter_generate_plot_mpl_from_workflow, sdec_regression_data
    ):
        _, plotter = plotter_generate_plot_mpl_from_workflow
        param_idx = plotter._param_idx
        regression_file = sdec_regression_data / f"test_generate_plot_mpl__plotter_generate_plot_ply{param_idx}__.h5"
        
        # Compare species names and color lists
        expected_species = pd.read_hdf(regression_file, key="plot_data_hdf/_species_name", mode='r')
        expected_colors = pd.read_hdf(regression_file, key="plot_data_hdf/_color_list", mode='r')
        
        np.testing.assert_array_equal(plotter._species_name, expected_species.values.flatten())
        
        color_list = [item for subitem in plotter._color_list for item in subitem]
        np.testing.assert_array_equal(color_list, expected_colors.values.flatten())

    def test_workflow_simulation_data_identical(self, plotter, plotter_from_workflow):
        # Calculate plotting data with identical parameters
        params = {
            "packets_mode": "virtual",
            "packet_wvl_range": [500, 9000] * u.AA,
            "distance": None,
            "nelements": None,
        }
        
        plotter._parse_species_list(None)
        plotter._calculate_plotting_data(**params)
        
        plotter_from_workflow._parse_species_list(None)
        plotter_from_workflow._calculate_plotting_data(**params)
        
        # Test emission luminosities are identical
        pd.testing.assert_frame_equal(
            plotter.emission_luminosities_df,
            plotter_from_workflow.emission_luminosities_df,
            atol=0, rtol=RELATIVE_TOLERANCE_SDEC
        )
        
        # Test absorption luminosities are identical  
        pd.testing.assert_frame_equal(
            plotter.absorption_luminosities_df,
            plotter_from_workflow.absorption_luminosities_df,
            atol=0, rtol=RELATIVE_TOLERANCE_SDEC
        )

    def test_mpl_image_workflow(self, plotter_generate_plot_mpl_from_workflow, tmp_path, sdec_regression_data):
        fig, plotter = plotter_generate_plot_mpl_from_workflow
        param_idx = plotter._param_idx
        
        # Save actual image
        actual_image_path = tmp_path / f"test_mpl_image_workflow_{param_idx}.png"
        fig.figure.savefig(actual_image_path)
        
        # Path to expected image
        expected_image_path = sdec_regression_data / f"test_mpl_image__plotter_generate_plot_mpl{param_idx}__.png"
        
        # Compare images
        compare_images(str(expected_image_path), str(actual_image_path), tol=1e-3)

    @pytest.mark.parametrize(
        "attribute", ["_full_species_list", "_species_list", "_keep_colour"]
    )
    def test_parse_species_list_workflow(self, plotter_from_workflow, attribute, sdec_regression_data):
        # Parse species list on workflow plotter
        plotter_from_workflow._parse_species_list(self.species_list[0])
        
        # Load expected data from regression files
        expected_file = sdec_regression_data / f"test_parse_species_list__{attribute}__.npy"
        expected_data = np.load(expected_file)
        
        actual_data = getattr(plotter_from_workflow, attribute)
        
        if attribute == "_full_species_list":
            np.testing.assert_equal(actual_data, expected_data)
        else:
            np.testing.assert_allclose(actual_data, expected_data, atol=0, rtol=RELATIVE_TOLERANCE_SDEC)
