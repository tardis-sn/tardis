from copy import deepcopy
from itertools import product

import astropy.units as u
import numpy as np
import pytest
from matplotlib.testing.compare import compare_images

from tardis.base import run_tardis
from tardis.io.util import HDFWriterMixin
from tardis.visualization.tools.liv_plot import LIVPlotter
from tardis.tests.fixtures.regression_data import RegressionData


class PlotDataHDF(HDFWriterMixin):
    def __init__(self, **kwargs):
        self.hdf_properties = []
        for key, value in kwargs.items():
            setattr(self, key, value)
            self.hdf_properties.append(key)


@pytest.fixture(scope="module")
def simulation_simple(config_verysimple, atomic_dataset):
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
        plot_data = PlotDataHDF(**property_group)
        return plot_data

    def test_generate_plot_mpl(
        self, generate_plot_mpl_hdf, plotter_generate_plot_mpl, request
    ):
        fig, _ = plotter_generate_plot_mpl
        regression_data = RegressionData(request)
        expected = regression_data.sync_hdf_store(generate_plot_mpl_hdf)
        for item in ["_species_name", "_color_list", "step_x", "step_y"]:
            np.testing.assert_array_equal(
                expected.get("plot_data_hdf/" + item).values.flatten(),
                getattr(generate_plot_mpl_hdf, item),
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

    @pytest.fixture(scope="function", params=combinations)
    def plotter_generate_plot_ply(self, request, plotter):
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
    def generate_plot_plotly_hdf(self, plotter_generate_plot_ply, request):
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
        plot_data = PlotDataHDF(**property_group)
        return plot_data

    def test_generate_plot_ply(
        self, generate_plot_plotly_hdf, plotter_generate_plot_ply, request
    ):
        fig, _ = plotter_generate_plot_ply
        regression_data = RegressionData(request)
        expected = regression_data.sync_hdf_store(generate_plot_plotly_hdf)

        for item in ["_species_name", "_color_list", "step_x", "step_y"]:
            np.testing.assert_array_equal(
                expected.get("plot_data_hdf/" + item).values.flatten(),
                getattr(generate_plot_plotly_hdf, item),
            )
