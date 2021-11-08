"""Tests for SDEC Plots."""
from tardis.base import run_tardis
from tardis.io.config_reader import Configuration
import pytest
import pandas as pd
import numpy as np
from copy import deepcopy
from tardis.visualization.tools.sdec_plot import SDECData, SDECPlotter
import astropy.units as u
import h5py


@pytest.fixture(scope="module")
def config_verysimple(atomic_dataset):
    config = Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
    )
    return config


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

    atomic_data = deepcopy(atomic_dataset)  # TODO: why deepcopy?
    sim = run_tardis(
        config_verysimple,
        atom_data=atomic_data,
        show_convergence_plots=False,
    )
    return sim


class TestSDECPlotter:
    @pytest.fixture(scope="class")
    def plotter(self, simulation_simple):
        return SDECPlotter.from_simulation(simulation_simple)

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, 50 * u.Mpc])
    @pytest.mark.parametrize("nelements", [1, 3])
    @pytest.mark.parametrize("species", [["Si II", "Ca II", "C", "Fe I-V"]])
    def test_calculate_plotting_data(
        self,
        request,
        simulation_simple,
        plotter,
        packets_mode,
        packet_wvl_range,
        distance,
        nelements,
        species,
    ):

        # each group is a different combination of arguments
        subgroup_name = request.node.callspec.id
        ref_data_path = ""
        # TODO: get refdata path

        if request.config.getoption("--generate-reference"):
            plotter._parse_species_list(species)
            plotter._calculate_plotting_data(
                packets_mode, packet_wvl_range, distance, nelements
            )
            # TODO: delete object at ref_data_path

            with h5py.File("sdec_ref.h5", "a") as file:
                group = file.create_group(subgroup_name)
                group.attrs["packet_wvl_range"] = packet_wvl_range.cgs.value
                group.attrs["packets_mode"] = packets_mode
                group.attrs["distance"] = distance.cgs.value
                group.attrs["nelements"] = nelements
                group.attrs["species"] = species

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
                pytest.skip(
                    f"SDEC test data saved at: {ref_data_path}",
                    allow_module_level=True,
                )
        else:
            # use the subgroup id to iterate over the hdf file
            plotter._parse_species_list(species)
            plotter._calculate_plotting_data(
                packets_mode, packet_wvl_range, distance, nelements=nelements
            )

            with h5py.File("sdec_ref.h5", "r") as file:
                group = file[subgroup_name]

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

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, 50 * u.Mpc])
    @pytest.mark.parametrize("show_modeled_spectrum", [True, False])
    def test_generate_plot_mpl(
        self,
        plotter,
        packets_mode,
        packet_wvl_range,
        distance,
        show_modeled_spectrum,
    ):
        plotter.generate_plot_mpl(
            packets_mode, packet_wvl_range, distance, show_modeled_spectrum
        )

    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, 50 * u.Mpc])
    @pytest.mark.parametrize("show_modeled_spectrum", [True, False])
    def test_generate_plot_ply(
        self,
        plotter,
        packets_mode,
        packet_wvl_range,
        distance,
        show_modeled_spectrum,
    ):
        fig = plotter.generate_plot_ply(
            packets_mode, packet_wvl_range, distance, show_modeled_spectrum
        )

        for trace in fig.data:
            if trace.name == f"{packets_mode.capitalize()} Spectrum":
                assert (trace.x == plotter.plot_wavelength.value).all()
                assert (
                    trace.y == plotter.modeled_spectrum_luminosity.value
                ).all()

            if trace.name == "Blackbody Photosphere":
                assert (trace.x == plotter.plot_wavelength.value).all()
                assert (trace.y == plotter.photosphere_luminosity.value).all()

            if trace.name == "No interaction":
                assert (
                    trace.x == plotter.emission_luminosities_df.index.values
                ).all()
                assert (
                    trace.y == plotter.emission_luminosities_df.noint.values
                ).all()

            if trace.name == "Electron Scatter Only":
                assert (
                    trace.x == plotter.emission_luminosities_df.index.values
                ).all()
                assert (
                    trace.y == plotter.emission_luminosities_df.escatter.values
                ).all()

            if trace.name == "Other Elements":
                assert (
                    trace.x == plotter.emission_luminosities_df.index.values
                ).all()
                assert (
                    trace.y == plotter.emission_luminosities_df.other.values
                ).all()

    # TODO: Call generate_plot_mpl() with several PnCs of parameters and test
    # almost plotter's properties (esp. those saved by calculate_plotting_data)
    # against saved test data for those combinations (need to figure out a good
    # structure for saving test data for different PnCs)
