"""Tests for SDEC Plots."""
from tardis.base import run_tardis
from tardis.io.config_reader import Configuration
import pytest
import pandas as pd
import numpy as np
from copy import deepcopy
from tardis.visualization.tools.sdec_plot import SDECData, SDECPlotter
import astropy.units as u

@pytest.fixture(scope="module")
def config_verysimple(atomic_dataset):
    # TODO: do we really need this?
    config = Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
    )
    return config


class TestSDECPlotter:
    @pytest.fixture(scope="class")
    def plotter(self, config_verysimple, atomic_dataset):
        """Instantiate SDEC plotter using a simple simulation model."""
        # Setup simulation configuration using config_verysimple and
        # override properties in such a way to make the simulation run faster
        config_verysimple.montecarlo.iterations = 3
        config_verysimple.montecarlo.no_of_packets = 4000
        config_verysimple.montecarlo.last_no_of_packets = -1
        config_verysimple.spectrum.virtual.virtual_packet_logging = True
        config_verysimple.spectrum.num = 2000

        atomic_data = deepcopy(atomic_dataset) # TODO: why deepcopy?
        sim = run_tardis(config_verysimple, atom_data=atomic_data, show_convergence_plots=False)
        return SDECPlotter.from_simulation(sim)
    
    @pytest.mark.parametrize("species_list", [["Si II", "Si I-V", "Ca"], None])
    def test_parse_species_list(self, plotter, species_list):
        """Test the parse species list function"""
        plotter._parse_species_list(species_list)
    
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, 50 * u.Mpc])
    @pytest.mark.parametrize("show_modeled_spectrum", [True, False])
    def test_generate_plot_mpl(self, plotter, packets_mode, packet_wvl_range, distance, show_modeled_spectrum):
        plotter.generate_plot_mpl(packets_mode, packet_wvl_range, distance, show_modeled_spectrum)
       
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, 50 * u.Mpc])
    @pytest.mark.parametrize("show_modeled_spectrum", [True, False]) 
    def test_generate_plot_ply(self, plotter, packets_mode, packet_wvl_range, distance, show_modeled_spectrum):
        fig = plotter.generate_plot_ply(packets_mode, packet_wvl_range, distance, show_modeled_spectrum)
        
        for trace in fig.data:
            if trace.name == f"{packets_mode.capitalize()} Spectrum":
                assert (trace.x == plotter.plot_wavelength.value).all()
                assert (trace.y == plotter.modeled_spectrum_luminosity.value).all()
            
            if trace.name == "Blackbody Photosphere":
                assert (trace.x == plotter.plot_wavelength.value).all()
                assert (trace.y == plotter.photosphere_luminosity.value).all()
            
            if trace.name == "No interaction":
                assert (trace.x == plotter.emission_luminosities_df.index.values).all()
                assert (trace.y == plotter.emission_luminosities_df.noint.values).all()
            
            if trace.name == "Electron Scatter Only":
                assert (trace.x == plotter.emission_luminosities_df.index.values).all()
                assert (trace.y == plotter.emission_luminosities_df.escatter.values).all()
                                
            
        
        
    # TODO: Call generate_plot_mpl() with several PnCs of parameters and test
    # almost plotter's properties (esp. those saved by calculate_plotting_data)
    # against saved test data for those combinations (need to figure out a good
    # structure for saving test data for different PnCs)
