"""Tests for LIV Plots."""

import os
from copy import deepcopy
import json

import numpy as np
import pytest
import tables
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

from tardis.base import run_tardis
from tardis.visualization.tools.liv_plot import LIVPlotter


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


@pytest.fixture(scope="module")
def liv_ref_data_path(tardis_ref_path):
    """
    Return the path to the reference data for the LIV plots.

    Parameters
    ----------
    tardis_ref_path : str
        Path to the reference data directory.

    Returns
    -------
    str
        Path to LIV reference data.
    """
    return os.path.abspath(os.path.join(tardis_ref_path, "liv_ref.h5"))


class TestLIVPlotter:
    """Test the LIVPlotter class."""

    @pytest.fixture(scope="class", autouse=True)
    def create_hdf_file(self, request, liv_ref_data_path):
        """
        Create an HDF5 file object.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        liv_ref_data_path : str
            Path to the reference data for the LIV plots.

        Yields
        -------
        h5py._hl.files.File
            HDF5 file object.
        """
        cls = type(self)
        if request.config.getoption("--generate-reference"):
            cls.hdf_file = tables.open_file(liv_ref_data_path, "w")
        else:
            cls.hdf_file = tables.open_file(liv_ref_data_path, "r")
        yield cls.hdf_file
        cls.hdf_file.close()

    @pytest.fixture(scope="class")
    def plotter(self, simulation_simple):
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

    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"]]
    )
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("nelements", [1, None])
    def test_parse_species_list(
        self, request, plotter, species_list, packets_mode, nelements
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
        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group, name="_species_list", obj=plotter._species_list
            )
            self.hdf_file.create_carray(
                group, name="_keep_colour", obj=plotter._keep_colour
            )
            species_mapped_json = json.dumps(
                convert_to_native_type(plotter._species_mapped)
            )
            self.hdf_file.create_array(
                group,
                name="_species_mapped",
                obj=np.array([species_mapped_json], dtype="S"),
            )

            pytest.skip("Reference data was generated during this run.")
        else:
            group = self.hdf_file.get_node("/" + subgroup_name)

            np.testing.assert_allclose(
                np.asarray(plotter._species_list),
                self.hdf_file.get_node(group, "_species_list"),
            )
            np.testing.assert_allclose(
                np.asarray(plotter._keep_colour),
                self.hdf_file.get_node(group, "_keep_colour"),
            )
            species_mapped_array = self.hdf_file.get_node(
                group, "_species_mapped"
            ).read()
            species_mapped_json = (
                species_mapped_array[0].decode()
                if isinstance(species_mapped_array[0], bytes)
                else species_mapped_array[0]
            )
            species_mapped_dict = json.loads(species_mapped_json)
            species_mapped_dict = {
                int(key): value for key, value in species_mapped_dict.items()
            }
            assert plotter._species_mapped == species_mapped_dict

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
    ):
        """
        Test _parse_species_list method.

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
        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group, name="plot_data", obj=plot_data_list
            )

            self.hdf_file.create_carray(
                group, name="plot_colors", obj=plotter.plot_colors
            )
            self.hdf_file.create_carray(
                group, name="new_bin_edges", obj=plotter.new_bin_edges
            )
            pytest.skip("Reference data was generated during this run.")

        else:
            group = self.hdf_file.get_node("/" + subgroup_name)

            np.testing.assert_allclose(
                np.asarray(plot_data_list),
                self.hdf_file.get_node(group, "plot_data"),
            )

            np.testing.assert_allclose(
                np.asarray(plotter.plot_colors),
                self.hdf_file.get_node(group, "plot_colors"),
            )
            np.testing.assert_allclose(
                np.asarray(plotter.new_bin_edges),
                self.hdf_file.get_node(group, "new_bin_edges"),
            )

    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"], None]
    )
    @pytest.mark.parametrize("nelements", [1, 2, 3, 4])
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("xlog_scale", [True, False])
    @pytest.mark.parametrize("ylog_scale", [True, False])
    @pytest.mark.parametrize("num_bins", [10, 25])
    @pytest.mark.parametrize("velocity_range", [(12500, 15000), (15050, 19000)])
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

        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group, name="_species_name", obj=plotter._species_name
            )
            self.hdf_file.create_carray(
                group, name="_color_list", obj=plotter._color_list
            )
            self.hdf_file.create_carray(
                group, name="step_x", obj=plotter.step_x
            )
            self.hdf_file.create_carray(
                group, name="step_y", obj=plotter.step_y
            )
            fig_subgroup = self.hdf_file.create_group(
                group,
                name="fig_data",
            )

            for index, data in enumerate(fig.get_children()):
                trace_group = self.hdf_file.create_group(
                    fig_subgroup,
                    name="_" + str(index),
                )
                if isinstance(data.get_label(), str):
                    self.hdf_file.create_array(
                        trace_group, name="label", obj=data.get_label().encode()
                    )

                # save artists which correspond to element contributions
                if isinstance(data, PolyCollection):
                    for index, path in enumerate(data.get_paths()):
                        self.hdf_file.create_carray(
                            trace_group,
                            name="path" + str(index),
                            obj=path.vertices,
                        )
                # save line plots
                if isinstance(data, Line2D):
                    self.hdf_file.create_carray(
                        trace_group,
                        name="data",
                        obj=data.get_xydata(),
                    )
                    self.hdf_file.create_carray(
                        trace_group, name="path", obj=data.get_path().vertices
                    )
            pytest.skip("Reference data was generated during this run.")

        else:
            group = self.hdf_file.get_node("/" + subgroup_name)

            assert (
                plotter._species_name
                == self.hdf_file.get_node(group, "_species_name")
                .read()
                .astype(str),
            )
            np.testing.assert_allclose(
                np.asarray(np.asarray(plotter._color_list)),
                self.hdf_file.get_node(group, "_color_list"),
            )
            np.testing.assert_allclose(
                np.asarray(plotter.step_x),
                self.hdf_file.get_node(group, "step_x"),
            )
            np.testing.assert_allclose(
                np.asarray(plotter.step_y),
                self.hdf_file.get_node(group, "step_y"),
            )
            fig_subgroup = self.hdf_file.get_node(group, "fig_data")
            for index, data in enumerate(fig.get_children()):
                trace_group = self.hdf_file.get_node(
                    fig_subgroup, "_" + str(index)
                )
                if isinstance(data.get_label(), str):
                    assert (
                        data.get_label()
                        == self.hdf_file.get_node(trace_group, "label")
                        .read()
                        .decode()
                    )

                # test element contributions
                if isinstance(data, PolyCollection):
                    for index, path in enumerate(data.get_paths()):
                        np.testing.assert_allclose(
                            path.vertices,
                            self.hdf_file.get_node(
                                trace_group, "path" + str(index)
                            ),
                        )
                # compare line plot data
                if isinstance(data, Line2D):
                    np.testing.assert_allclose(
                        data.get_xydata(),
                        self.hdf_file.get_node(trace_group, "data"),
                    )
                    np.testing.assert_allclose(
                        data.get_path().vertices,
                        self.hdf_file.get_node(trace_group, "path"),
                    )

    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"], None]
    )
    @pytest.mark.parametrize("nelements", [1, 2, 3, 4])
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("xlog_scale", [True, False])
    @pytest.mark.parametrize("ylog_scale", [True, False])
    @pytest.mark.parametrize("num_bins", [10, 25])
    @pytest.mark.parametrize("velocity_range", [(12500, 15000), (15050, 25000)])
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

        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group, name="_species_name", obj=plotter._species_name
            )
            self.hdf_file.create_carray(
                group, name="_color_list", obj=plotter._color_list
            )
            self.hdf_file.create_carray(
                group, name="step_x", obj=plotter.step_x
            )
            self.hdf_file.create_carray(
                group, name="step_y", obj=plotter.step_y
            )
            fig_subgroup = self.hdf_file.create_group(
                group,
                name="fig_data",
            )
            for index, data in enumerate(fig.data):
                trace_group = self.hdf_file.create_group(
                    fig_subgroup,
                    name="_" + str(index),
                )
                if data.stackgroup:
                    self.hdf_file.create_array(
                        trace_group,
                        name="stackgroup",
                        obj=data.stackgroup.encode(),
                    )
                if data.name:
                    self.hdf_file.create_array(
                        trace_group,
                        name="name",
                        obj=data.name.encode(),
                    )
                self.hdf_file.create_carray(
                    trace_group,
                    name="x",
                    obj=data.x,
                )
                self.hdf_file.create_carray(
                    trace_group,
                    name="y",
                    obj=data.y,
                )
            pytest.skip("Reference data was generated during this run.")

        else:
            group = self.hdf_file.get_node("/", subgroup_name)

            assert (
                plotter._species_name
                == self.hdf_file.get_node(group, "_species_name")
                .read()
                .astype(str),
            )
            # test output of the _make_colorbar_colors function
            np.testing.assert_allclose(
                np.asarray(np.asarray(plotter._color_list)),
                self.hdf_file.get_node(group, "_color_list"),
            )
            np.testing.assert_allclose(
                np.asarray(plotter.step_x),
                self.hdf_file.get_node(group, "step_x"),
            )
            np.testing.assert_allclose(
                np.asarray(plotter.step_y),
                self.hdf_file.get_node(group, "step_y"),
            )
            fig_subgroup = self.hdf_file.get_node(group, "fig_data")
            for index, data in enumerate(fig.data):
                trace_group = self.hdf_file.get_node(
                    fig_subgroup, "_" + str(index)
                )
                if data.stackgroup:
                    assert (
                        data.stackgroup
                        == self.hdf_file.get_node(trace_group, "stackgroup")
                        .read()
                        .decode()
                    )
                if data.name:
                    assert (
                        data.name
                        == self.hdf_file.get_node(trace_group, "name")
                        .read()
                        .decode()
                    )
                np.testing.assert_allclose(
                    self.hdf_file.get_node(trace_group, "x"), data.x
                )
                np.testing.assert_allclose(
                    self.hdf_file.get_node(trace_group, "y"), data.y
                )
