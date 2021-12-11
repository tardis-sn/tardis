"""Tests for SDEC Plots."""
from tardis.base import run_tardis
import pytest
import pandas as pd
import numpy as np
import os
from copy import deepcopy
from tardis.visualization.tools.sdec_plot import SDECData, SDECPlotter
import astropy.units as u
import h5py
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D


@pytest.fixture(scope="module")
def simulation_simple(config_verysimple, atomic_dataset):
    """Instantiate SDEC plotter using a simple simulation model."""
    # Setup simulation configuration using config_verysimple and
    # override properties in such a way to make the simulation run faster
    config_verysimple.montecarlo.iterations = 3
    config_verysimple.montecarlo.no_of_packets = 4000
    config_verysimple.montecarlo.last_no_of_packets = -1
    config_verysimple.spectrum.virtual.virtual_packet_logging = True
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
    return os.path.abspath(os.path.join(tardis_ref_path, "sdec_ref.h5"))

class TestSDECPlotter(object):
    @classmethod
    @pytest.fixture(scope="class", autouse=True)
    def create_hdf_file(self, request, sdec_ref_data_path):
        cls = type(self)
        if request.config.getoption("--generate-reference"):
            cls.hdf_file = h5py.File(sdec_ref_data_path, "w")
        else:
            cls.hdf_file = h5py.File(sdec_ref_data_path, "r")
        yield cls.hdf_file
        cls.hdf_file.close()
    
    @pytest.fixture(scope="class")
    def plotter(self, simulation_simple):
        return SDECPlotter.from_simulation(simulation_simple)

    @pytest.fixture(scope="class")
    def observed_spectrum(self):
        test_data = np.loadtxt(
            "tardis/visualization/tools/tests/data/observed_spectrum_test_data.dat"
        )
        observed_spectrum_wavelength, observed_spectrum_flux = test_data.T
        observed_spectrum_wavelength = observed_spectrum_wavelength * u.AA
        observed_spectrum_flux = (
            observed_spectrum_flux * u.erg / (u.s * u.cm ** 2 * u.AA)
        )
        return observed_spectrum_wavelength, observed_spectrum_flux

    @pytest.mark.parametrize("species", [["Si II", "Ca II", "C", "Fe I-V"]])
    def test_parse_species_list(self, request, plotter, species):
        plotter._parse_species_list(species)
        subgroup_name = request.node.callspec.id
        
        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(subgroup_name)
            group.create_dataset(
                "_full_species_list", data=plotter._full_species_list
            )
            group.create_dataset("_species_list", data=plotter._species_list)
            group.create_dataset("_keep_colour", data=plotter._keep_colour)
            # pytest skip this test?
        else:
            group = self.hdf_file[subgroup_name]
            assert (
                plotter._full_species_list
                == group.get("_full_species_list").asstr()[()].tolist()
            )
            np.testing.assert_allclose(
                np.asarray(plotter._species_list),
                group.get("_species_list")[()],
            )
            np.testing.assert_allclose(
                np.asarray(plotter._keep_colour), group.get("_keep_colour")[()]
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
        plotter._calculate_plotting_data(
            packets_mode, packet_wvl_range, distance, nelements
        )

        # each group is a different combination of arguments
        subgroup_name = request.node.callspec.id

        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(subgroup_name)
            group.attrs["packet_wvl_range"] = packet_wvl_range.cgs.value
            group.attrs["packets_mode"] = packets_mode
            group.attrs["distance"] = distance.cgs.value
            group.attrs["nelements"] = nelements

            group.create_dataset(
                "plot_frequency_bins",
                data=plotter.plot_frequency_bins.cgs.value,
            )
            group.create_dataset(
                "plot_wavelength", data=plotter.plot_wavelength.cgs.value
            )
            group.create_dataset(
                "plot_frequency", data=plotter.plot_frequency.cgs.value
            )
            group.create_dataset(
                "packet_wvl_range_mask", data=plotter.packet_wvl_range_mask
            )
            group.create_dataset(
                "emission_species", data=plotter.emission_species
            )
            group.create_dataset(
                "absorption_species", data=plotter.absorption_species
            )
            group.create_dataset(
                "modeled_spectrum_luminosity",
                data=plotter.modeled_spectrum_luminosity.cgs.value,
            )
            group.create_dataset("lum_to_flux", data=plotter.lum_to_flux)
            group.create_dataset(
                "species", data=plotter.species.astype(np.float64)
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
                self.hdf_file.filename, key=f"{subgroup_name}/total_luminosities_df"
            )

            # pytest.skip(
            #     f"SDEC test data saved at: {ref_data_path}",
            #     allow_module_level=True,
            # )
        else:
            # use the subgroup id to iterate over the hdf file
            group = self.hdf_file[subgroup_name]

            np.testing.assert_allclose(
                plotter.plot_frequency_bins.cgs.value,
                group.get("plot_frequency_bins")[()],
            )
            np.testing.assert_allclose(
                plotter.plot_wavelength.cgs.value,
                group.get("plot_wavelength")[()],
            )
            np.testing.assert_allclose(
                plotter.plot_frequency.cgs.value,
                group.get("plot_frequency")[()],
            )
            np.testing.assert_allclose(
                plotter.modeled_spectrum_luminosity.cgs.value,
                group.get("modeled_spectrum_luminosity")[()],
            )

            np.testing.assert_allclose(
                plotter.packet_wvl_range_mask,
                group.get("packet_wvl_range_mask")[()],
            )
            np.testing.assert_allclose(
                plotter.absorption_species,
                group.get("absorption_species")[()],
            )
            np.testing.assert_allclose(
                plotter.emission_species, group.get("emission_species")[()]
            )
            np.testing.assert_allclose(
                plotter.species.astype(np.float64), group.get("species")[()]
            )

            if isinstance(plotter.lum_to_flux, u.quantity.Quantity):
                assert (
                    plotter.lum_to_flux.cgs.value
                    == group.get("lum_to_flux")[()]
                )
            else:
                assert plotter.lum_to_flux == group.get("lum_to_flux")[()]

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
    def test_generate_plot_mpl(
        self,
        request,
        plotter,
        packets_mode,
        packet_wvl_range,
        distance,
        show_modeled_spectrum,
        observed_spectrum,
    ):
        subgroup_name = "mpl" + request.node.callspec.id

        fig = plotter.generate_plot_mpl(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            show_modeled_spectrum=show_modeled_spectrum,
            observed_spectrum=observed_spectrum if distance else None,
        )

        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(subgroup_name)
            group.create_dataset(
                "_species_name", data=plotter._species_name
            )
            group.create_dataset("_color_list", data=plotter._color_list)

            fig_subgroup = group.create_group("fig_data")

            for index, data in enumerate(fig.get_children()):
                trace_group = fig_subgroup.create_group(str(index))

                if isinstance(data.get_label(), str):
                    trace_group.attrs["name"] = data.get_label()

                # save artists which correspond to element contributions
                if isinstance(data, PolyCollection):
                    for index, path in enumerate(data.get_paths()):
                        trace_group.create_dataset(
                            "path" + str(index),
                            data=path.vertices,
                        )
                # save line plots
                if isinstance(data, Line2D):
                    trace_group.create_dataset(
                        "data", data=data.get_xydata()
                    )
                    trace_group.create_dataset(
                        "path", data=data.get_path().vertices
                    )

        else:
            group = self.hdf_file[subgroup_name]
            # test output of the _make_colorbar_labels function
            assert (
                plotter._species_name
                == group.get("_species_name").asstr()[()].tolist()
            )
            # test output of the _make_colorbar_colors function
            np.testing.assert_allclose(
                np.asarray(np.asarray(plotter._color_list)),
                group.get("_color_list")[()],
            )
            fig_subgroup = group.get("fig_data")
            for index, data in enumerate(fig.get_children()):
                trace_group = fig_subgroup.get(str(index))

                if isinstance(data.get_label(), str):
                    assert trace_group.attrs["name"] == data.get_label()

                # test element contributions
                if isinstance(data, PolyCollection):
                    for index, path in enumerate(data.get_paths()):
                        np.testing.assert_allclose(
                            trace_group.get("path" + str(index))[()],
                            path.vertices,
                        )
                # compare line plot data
                if isinstance(data, Line2D):
                    np.testing.assert_allclose(
                        trace_group.get("data")[()], data.get_xydata()
                    )
                    np.testing.assert_allclose(
                        trace_group.get("path")[()],
                        data.get_path().vertices,
                    )

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, 50 * u.Mpc])
    @pytest.mark.parametrize("show_modeled_spectrum", [True, False])
    def test_generate_plot_ply(
        self,
        request,
        plotter,
        packets_mode,
        packet_wvl_range,
        distance,
        show_modeled_spectrum,
        observed_spectrum
    ):
        subgroup_name = "ply" + request.node.callspec.id

        fig = plotter.generate_plot_ply(
            packets_mode=packets_mode,
            packet_wvl_range=packet_wvl_range,
            distance=distance,
            show_modeled_spectrum=show_modeled_spectrum,
            observed_spectrum=observed_spectrum if distance else None,
        )

        if request.config.getoption("--generate-reference"):
            group = self.hdf_file.create_group(subgroup_name)
            group.create_dataset(
                "_species_name", data=plotter._species_name
            )
            group.create_dataset("_color_list", data=plotter._color_list)

            fig_subgroup = group.create_group("fig_data")
            for index, data in enumerate(fig.data):
                trace_group = fig_subgroup.create_group(str(index))
                if data.name:
                    trace_group.attrs["name"] = data.name
                if data.stackgroup:
                    trace_group.attrs["stackgroup"] = data.stackgroup
                trace_group.create_dataset("x", data=data.x)
                trace_group.create_dataset("y", data=data.y)
        else:
            group = self.hdf_file[subgroup_name]
            # test output of the _make_colorbar_labels function
            assert (
                plotter._species_name
                == group.get("_species_name").asstr()[()].tolist()
            )
            # test output of the _make_colorbar_colors function
            np.testing.assert_allclose(
                np.asarray(np.asarray(plotter._color_list)),
                group.get("_color_list")[()],
            )

            fig_subgroup = group.get("fig_data")
            for index, data in enumerate(fig.data):
                trace_group = fig_subgroup.get(str(index))
                if data.name:
                    assert data.name == trace_group.attrs["name"]
                if data.stackgroup:
                    assert data.stackgroup == trace_group.attrs["stackgroup"]

                np.testing.assert_allclose(trace_group.get("x")[()], data.x)
                np.testing.assert_allclose(trace_group.get("y")[()], data.y)

    # TODO: Call generate_plot_mpl() with several PnCs of parameters and test
    # almost plotter's properties (esp. those saved by calculate_plotting_data)
    # against saved test data for those combinations (need to figure out a good
    # structure for saving test data for different PnCs)
