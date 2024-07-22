"""Tests for LIV Plots."""

import os
from copy import deepcopy

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
import tables
from matplotlib.collections import PolyCollection

from tardis.base import run_tardis
from tardis.visualization.tools.liv_plot import LIVPlotter
from matplotlib.collections import PolyCollection


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
        show_convergence_plots=False,
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

    @pytest.mark.parametrize("species", [["Si II", "Ca II", "C", "Fe I-V"]])
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("nelements", [1, None])
    def test_parse_species_list(
        self, request, plotter, species, packets_mode, nelements
    ):
        """
        Test _parse_species_list method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.liv_plot.LIVPlotter
        species : list
        """
        plotter._parse_species_list(species)
        subgroup_name = make_valid_name(request.node.callspec.id)
        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group, name="_full_species_list", obj=plotter._full_species_list
            )
            self.hdf_file.create_carray(
                group, name="_species_list", obj=plotter._species_list
            )
            self.hdf_file.create_carray(
                group, name="_keep_colour", obj=plotter._keep_colour
            )
            self.hdf_file.create_carray(
                group, name="species_mapped", obj=plotter.species_mapped
            )
            pytest.skip("Reference data was generated during this run.")
        else:
            group = self.hdf_file.get_node("/" + subgroup_name)

            np.testing.assert_equal(
                np.asarray(plotter._full_species_list),
                self.hdf_file.get_node(group, "_full_species_list")
                .read()
                .astype(str),
            )

            np.testing.assert_allclose(
                np.asarray(plotter._species_list),
                self.hdf_file.get_node(group, "_species_list"),
            )
            np.testing.assert_allclose(
                np.asarray(plotter._keep_colour),
                self.hdf_file.get_node(group, "_keep_colour"),
            )
            np.testing.assert_equal(
                np.asarray(plotter.species_mapped),
                self.hdf_file.get_node(group, "species_mapped").read(),
            )

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    def test_generate_plot_data(self, request, plotter, packets_mode):
        """
        Test generate_plot_data method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.liv_plot.LIVPlotter
        packets_mode : str
        """
        plotter._generate_plot_data(packets_mode)

        subgroup_name = make_valid_name(request.node.callspec.id)
        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group, name="plot_data", obj=plotter.plot_data
            )

            self.hdf_file.create_carray(
                group, name="plot_color", obj=plotter.plot_color
            )
            pytest.skip("Reference data was generated during this run.")
        else:
            group = self.hdf_file.get_node("/" + subgroup_name)

            np.testing.assert_allclose(
                np.asarray(plotter.plot_data),
                self.hdf_file.get_node(group, "plot_data"),
            )

            np.testing.assert_allclose(
                np.asarray(plotter.plot_color),
                self.hdf_file.get_node(group, "plot_color"),
            )

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"]]
    )
    @pytest.mark.parametrize("num_bins", [5, 10, 25, 40])
    @pytest.mark.parametrize("nelements", [1, None])
    def test_prepare_plot_data(
        self,
        request,
        plotter,
        packets_mode,
        species_list,
        num_bins,
        nelements,
    ):
        plotter._prepare_plot_data(
            packets_mode, species_list, num_bins, nelements
        )

        subgroup_name = make_valid_name(request.nod.callspec.id)
        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )

            self.hdf_file.create_carray(
                group, name="new_bin_edges", obj=plotter.new_bin_edges
            )
            pytest.skip("Reference data was generated during this run.")

        else:
            group = self.hdf_file.get_node("/" + subgroup_name)

            np.testing.assert_allclose(
                np.asarray(plotter.new_bin_edges),
                self.hdf_file.get_node(group, "new_bin_edges"),
            )

    def test_get_step_plot_data(self, request, plotter):
        """
        Test get_step_plot_data method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.liv_plot.LIVPlotter
        """
        step_plot_data = plotter.get_step_plot_data()

        subgroup_name = make_valid_name(request.node.callspec.id)
        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group, name="step_x", obj=plotter.step_x
            )
            self.hdf_file.create_carray(
                group, name="step_y", obj=plotter.step_y
            )
            pytest.skip("Reference data was generated during this run.")
        else:
            group = self.hdf_file.get_node("/" + subgroup_name)
            np.testing.assert_allclose(
                np.asarray(plotter.step_x),
                self.hdf_file.get_node(group, "step_x"),
            )
            np.testing.assert_allclose(
                np.asarray(plotter.step_y),
                self.hdf_file.get_node(group, "step_y"),
            )

    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"]]
    )
    @pytest.mark.parametrize("nelements", [1, None])
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("xlog_scale", [True, False])
    @pytest.mark.parametrize("ylog_scale", [True, False])
    @pytest.mark.parametrize("num_bins", [5, 10, 25, 40])
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
        """
        subgroup_name = make_valid_name("mpl" + request.node.callspec.id)
        fig = plotter.generate_plot_mpl(
            species_list,
            nelements,
            packets_mode,
            xlog_scale,
            ylog_scale,
            num_bins,
            velocity_range,
        )

        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group,
                name="mpl_fig",
                obj=fig,
            )
            pytest.skip("Reference data was generated during this run.")
        else:
            group = self.hdf_file.get_node("/" + subgroup_name)
            np.testing.assert_allclose(
                fig,
                self.hdf_file.get_node(group, "mpl_fig"),
            )

    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"]]
    )
    @pytest.mark.parametrize("nelements", [1, 2, 3, 4])
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("xlog_scale", [True, False])
    @pytest.mark.parametrize("ylog_scale", [True, False])
    @pytest.mark.parametrize("num_bins", [5, 10, 25, 40])
    @pytest.mark.parametrize("velocity_range", [(12500, 15000), (15050, 19000)])
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
        """
        subgroup_name = make_valid_name("ply" + request.node.callspec.id)
        fig = plotter.generate_plot_ply(
            species_list,
            nelements,
            packets_mode,
            xlog_scale,
            ylog_scale,
            num_bins,
            velocity_range,
        )

        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )
            self.hdf_file.create_carray(
                group,
                name="ply_fig",
                obj=fig,
            )
            pytest.skip("Reference data was generated during this run.")
        else:
            group = self.hdf_file.get_node("/" + subgroup_name)
            np.testing.assert_allclose(
                fig,
                self.hdf_file.get_node(group, "ply_fig"),
            )
