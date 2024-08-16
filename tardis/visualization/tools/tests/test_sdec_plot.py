"""Tests for SDEC Plots."""
from copy import deepcopy
from itertools import product

import astropy
import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from matplotlib.testing.compare import compare_images
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

from tardis.base import run_tardis
from tardis.io.util import HDFWriterMixin
from tardis.tests.fixtures.regression_data import RegressionData
from tardis.visualization.tools.sdec_plot import SDECPlotter


class PlotDataHDF(HDFWriterMixin):
    def __init__(self, **kwargs):
        self.hdf_properties = []
        for key, value in kwargs.items():
            setattr(self, key, value)
            self.hdf_properties.append(key)


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


@pytest.fixture(scope="module")
def simulation_simple(config_verysimple, atomic_dataset):
    """
    Instantiate SDEC plotter using a simple simulation model.

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
    # Setup simulation configuration using config_verysimple and
    # override properties in such a way to make the simulation run faster
    config_verysimple.montecarlo.iterations = 3
    config_verysimple.montecarlo.no_of_packets = 4000
    config_verysimple.montecarlo.last_no_of_packets = -1
    config_verysimple.spectrum.virtual.virtual_packet_logging = True
    config_verysimple.montecarlo.no_of_virtual_packets = 1
    config_verysimple.spectrum.num = 2000
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config_verysimple,
        atom_data=atomic_data,
        show_convergence_plots=False,
        log_level="CRITICAl",
    )
    return sim


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
    def plotter(self, simulation_simple, request):
        """
        Create a SDECPlotter object.

        Parameters
        ----------
        simulation_simple : tardis.simulation.base.Simulation
            Simulation object.

        Returns
        -------
        tardis.visualization.tools.sdec_plot.SDECPlotter
        """
        return SDECPlotter.from_simulation(simulation_simple)

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
    def test_parse_species_list(self, request, plotter, attribute):
        """
        Test _parse_species_list method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.sdec_plot.SDECPlotter
        species : list
        """
        # THIS NEEDS TO BE RUN FIRST. NOT INDEPENDENT TESTS
        plotter._parse_species_list(self.species_list[0])
        regression_data = RegressionData(request)
        data = regression_data.sync_ndarray(getattr(plotter, attribute))
        if attribute == "_full_species_list":
            np.testing.assert_equal(getattr(plotter, attribute), data)
        else:
            np.testing.assert_allclose(getattr(plotter, attribute), data)

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
        plotter._calculate_plotting_data(
            packets_mode, packet_wvl_range, distance, nelements
        )
        return plotter

    @pytest.fixture(scope="class")
    def calculate_plotting_data_hdf(
        self, request, plotter_calculate_plotting_data
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
        request,
    ):
        regression_data = RegressionData(request)
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
                    plot_object, expected.get(group + attribute_name)
                )
            if attribute_type == "attributes_pd":
                pd.testing.assert_frame_equal(
                    plot_object, expected.get(group + attribute_name)
                )

    @pytest.fixture(scope="class", params=combinations)
    def plotter_generate_plot_mpl(self, request, observed_spectrum, plotter):
        (
            distance,
            packet_wvl_range,
            species_list,
            packets_mode,
            nelements,
            show_modeled_spectrum,
        ) = request.param
        if distance is None:
            observed_spectrum = None

        fig = plotter.generate_plot_mpl(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            show_modeled_spectrum=show_modeled_spectrum,
            observed_spectrum=observed_spectrum,
            nelements=nelements,
            species_list=species_list,
        )
        return fig, plotter

    @pytest.fixture(scope="class")
    def generate_plot_mpl_hdf(self, plotter_generate_plot_mpl, request):
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
        self, generate_plot_mpl_hdf, plotter_generate_plot_mpl, request
    ):
        fig, _ = plotter_generate_plot_mpl
        regression_data = RegressionData(request)
        expected = regression_data.sync_hdf_store(generate_plot_mpl_hdf)
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
                )
                np.testing.assert_allclose(
                    data.get_path().vertices,
                    expected.get("plot_data_hdf/" + "linepath" + str(index1)),
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

    @pytest.fixture(scope="class", params=combinations)
    def plotter_generate_plot_ply(self, request, observed_spectrum, plotter):
        (
            distance,
            packet_wvl_range,
            species_list,
            packets_mode,
            nelements,
            show_modeled_spectrum,
        ) = request.param
        if distance is None:
            observed_spectrum = None

        fig = plotter.generate_plot_ply(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            show_modeled_spectrum=show_modeled_spectrum,
            observed_spectrum=observed_spectrum,
            nelements=nelements,
            species_list=species_list,
        )
        return fig, plotter

    @pytest.fixture(scope="class")
    def generate_plot_plotly_hdf(self, plotter_generate_plot_ply, request):
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

    def test_generate_plot_mpl(
        self, generate_plot_plotly_hdf, plotter_generate_plot_ply, request
    ):
        fig, _ = plotter_generate_plot_ply
        regression_data = RegressionData(request)
        expected = regression_data.sync_hdf_store(generate_plot_plotly_hdf)

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
                data.x, expected.get(group + "x").values.flatten()
            )
            np.testing.assert_allclose(
                data.y, expected.get(group + "y").values.flatten()
            )

    def test_mpl_image(self, plotter_generate_plot_mpl, tmp_path, request):
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
