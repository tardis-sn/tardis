"""Tests for SDEC Plots."""
import os
from copy import deepcopy

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
import tables
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

from tardis.base import run_tardis
from tardis.visualization.tools.sdec_plot import SDECPlotter


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
    )
    return sim


@pytest.fixture(scope="module")
def sdec_ref_data_path(tardis_ref_path):
    """
    Return the path to the reference data for the SDEC plots.

    Parameters
    ----------
    tardis_ref_path : str
        Path to the reference data directory.

    Returns
    -------
    str
        Path to SDEC reference data.
    """
    return os.path.abspath(os.path.join(tardis_ref_path, "sdec_ref.h5"))


class TestSDECPlotter:
    """Test the SDECPlotter class."""

    @pytest.fixture(scope="class", autouse=True)
    def create_hdf_file(self, request, sdec_ref_data_path):
        """
        Create an HDF5 file object.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        sdec_ref_data_path : str
            Path to the reference data for the SDEC plots.

        Yields
        -------
        h5py._hl.files.File
            HDF5 file object.
        """
        cls = type(self)
        if request.config.getoption("--generate-reference"):
            cls.hdf_file = tables.open_file(sdec_ref_data_path, "w")

        else:
            cls.hdf_file = tables.open_file(sdec_ref_data_path, "r")
        yield cls.hdf_file
        cls.hdf_file.close()

    @pytest.fixture(scope="class")
    def plotter(self, simulation_simple):
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

    @pytest.mark.parametrize("species", [["Si II", "Ca II", "C", "Fe I-V"]])
    def test_parse_species_list(self, request, plotter, species):
        """
        Test _parse_species_list method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.sdec_plot.SDECPlotter
        species : list
        """
        # THIS NEEDS TO BE RUN FIRST. NOT INDEPENDENT TESTS
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
            pytest.skip("Reference data was generated during this run.")

        else:
            group = self.hdf_file.get_node("/" + subgroup_name)

            # because plotter._full_species_list is an array of strings
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

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, 50 * u.Mpc])
    @pytest.mark.parametrize("nelements", [1, 3])
    def test_calculate_plotting_data(
        self,
        request,
        plotter,
        packets_mode,
        packet_wvl_range,
        distance,
        nelements,
    ):
        """
        Test _calculate_plotting_data method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.sdec_plot.SDECPlotter
        packets_mode : str
        packet_wvl_range : astropy.units.quantity.Quantity
        distance : astropy.units.quantity.Quantity
        nelements : int
        """
        plotter._calculate_plotting_data(
            packets_mode, packet_wvl_range, distance, nelements
        )

        # each group is a different combination of arguments
        subgroup_name = make_valid_name(request.node.callspec.id)
        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(
                self.hdf_file.root,
                name=subgroup_name,
            )

            self.hdf_file.create_carray(
                group,
                name="plot_frequency_bins",
                obj=plotter.plot_frequency_bins.cgs.value,
            )

            self.hdf_file.create_carray(
                group,
                name="plot_wavelength",
                obj=plotter.plot_wavelength.cgs.value,
            )

            self.hdf_file.create_carray(
                group,
                name="plot_frequency",
                obj=plotter.plot_frequency.cgs.value,
            )

            self.hdf_file.create_carray(
                group,
                name="packet_wvl_range_mask",
                obj=plotter.packet_wvl_range_mask,
            )

            self.hdf_file.create_carray(
                group, name="emission_species", obj=plotter.emission_species
            )

            self.hdf_file.create_carray(
                group, name="absorption_species", obj=plotter.absorption_species
            )

            self.hdf_file.create_carray(
                group,
                name="modeled_spectrum_luminosity",
                obj=plotter.modeled_spectrum_luminosity.cgs.value,
            )

            if isinstance(plotter.lum_to_flux, u.quantity.Quantity):
                self.hdf_file.create_array(
                    group, name="lum_to_flux", obj=plotter.lum_to_flux.cgs.value
                )
            else:
                self.hdf_file.create_array(
                    group, name="lum_to_flux", obj=plotter.lum_to_flux
                )

            self.hdf_file.create_carray(
                group, name="species", obj=plotter.species.astype(np.float64)
            )

            plotter.absorption_luminosities_df.to_hdf(
                self.hdf_file.filename,
                key=f"{subgroup_name}/absorption_luminosities_df",
            )
            plotter.emission_luminosities_df.to_hdf(
                self.hdf_file.filename,
                key=f"{subgroup_name}/emission_luminosities_df",
            )
            plotter.total_luminosities_df.to_hdf(
                self.hdf_file.filename,
                key=f"{subgroup_name}/total_luminosities_df",
            )

            pytest.skip("Reference data was generated during this run.")

        else:
            # use the subgroup id to iterate over the hdf file
            group = self.hdf_file.get_node("/" + subgroup_name)

            np.testing.assert_allclose(
                plotter.plot_frequency_bins.cgs.value,
                self.hdf_file.get_node(group, "plot_frequency_bins"),
            )

            np.testing.assert_allclose(
                plotter.plot_wavelength.cgs.value,
                self.hdf_file.get_node(group, "plot_wavelength"),
            )

            np.testing.assert_allclose(
                plotter.plot_frequency.cgs.value,
                self.hdf_file.get_node(group, "plot_frequency"),
            )

            np.testing.assert_allclose(
                plotter.modeled_spectrum_luminosity.cgs.value,
                self.hdf_file.get_node(group, "modeled_spectrum_luminosity"),
            )

            np.testing.assert_allclose(
                plotter.packet_wvl_range_mask,
                self.hdf_file.get_node(group, "packet_wvl_range_mask"),
            )
            np.testing.assert_allclose(
                plotter.absorption_species,
                self.hdf_file.get_node(group, "absorption_species"),
            )

            np.testing.assert_allclose(
                plotter.emission_species,
                self.hdf_file.get_node(group, "emission_species"),
            )

            if isinstance(plotter.lum_to_flux, u.quantity.Quantity):
                assert (
                    plotter.lum_to_flux.cgs.value
                    == self.hdf_file.get_node(group, "lum_to_flux"),
                )
            else:
                assert plotter.lum_to_flux == self.hdf_file.get_node(
                    group, "lum_to_flux"
                )
            np.testing.assert_allclose(
                plotter.species.astype(np.float64),
                self.hdf_file.get_node(group, "species"),
            )
            pd.testing.assert_frame_equal(
                plotter.absorption_luminosities_df,
                pd.read_hdf(
                    self.hdf_file.filename,
                    key=f"{subgroup_name}/absorption_luminosities_df",
                ),
            )
            pd.testing.assert_frame_equal(
                plotter.emission_luminosities_df,
                pd.read_hdf(
                    self.hdf_file.filename,
                    key=f"{subgroup_name}/emission_luminosities_df",
                ),
            )
            pd.testing.assert_frame_equal(
                plotter.total_luminosities_df,
                pd.read_hdf(
                    self.hdf_file.filename,
                    key=f"{subgroup_name}/total_luminosities_df",
                ),
            )

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA, None])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, None])
    @pytest.mark.parametrize("show_modeled_spectrum", [True, False])
    @pytest.mark.parametrize("nelements", [1, None])
    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"], None]
    )
    def test_generate_plot_mpl(
        self,
        request,
        plotter,
        packets_mode,
        packet_wvl_range,
        distance,
        show_modeled_spectrum,
        observed_spectrum,
        nelements,
        species_list,
    ):
        """
        Test generate_plot_mpl method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.sdec_plot.SDECPlotter
        packets_mode : str
        packet_wvl_range : astropy.units.quantity.Quantity
        distance : astropy.units.quantity.Quantity
        show_modeled_spectrum : bool
        observed_spectrum : tuple of two astropy.units.quantity.Quantity values
        nelements : int
        species_list : list of str
        """
        subgroup_name = make_valid_name("mpl" + request.node.callspec.id)
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
            # test output of the _make_colorbar_labels function
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

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA, None])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, None])
    @pytest.mark.parametrize("show_modeled_spectrum", [True, False])
    @pytest.mark.parametrize("nelements", [1, None])
    @pytest.mark.parametrize(
        "species_list", [["Si II", "Ca II", "C", "Fe I-V"], None]
    )
    def test_generate_plot_ply(
        self,
        request,
        plotter,
        packets_mode,
        packet_wvl_range,
        distance,
        show_modeled_spectrum,
        observed_spectrum,
        nelements,
        species_list,
    ):
        """
        Test generate_plot_mpl method.

        Parameters
        ----------
        request : _pytest.fixtures.SubRequest
        plotter : tardis.visualization.tools.sdec_plot.SDECPlotter
        packets_mode : str
        packet_wvl_range : astropy.units.quantity.Quantity
        distance : astropy.units.quantity.Quantity
        show_modeled_spectrum : bool
        observed_spectrum : tuple of two astropy.units.quantity.Quantity values
        nelements : int
        species_list : list of str
        """
        subgroup_name = make_valid_name("ply" + request.node.callspec.id)
        if distance is not None:
            observed_spectrum = observed_spectrum
        else:
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
            # test output of the _make_colorbar_labels function
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
