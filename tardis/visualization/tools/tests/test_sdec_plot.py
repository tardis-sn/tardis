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
    def simulation_simple(self, config_verysimple, atomic_dataset):
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
        return sim
    
    @pytest.fixture(scope="class")
    def plotter(self, simulation_simple):
        return SDECPlotter.from_simulation(simulation_simple)
    
    
    
    @pytest.mark.parametrize("packets_mode", ["virtual", "real"])
    @pytest.mark.parametrize("packet_wvl_range", [[500, 9000] * u.AA, None])
    @pytest.mark.parametrize("distance", [10 * u.Mpc, 50 * u.Mpc])
    @pytest.mark.parametrize("nelements", [None])
    def test_calculate_plotting_data(self, simulation_simple,plotter, packets_mode, packet_wvl_range, distance, nelements):
        
        plotter._parse_species_list(None)
        plotter._calculate_plotting_data(packets_mode, packet_wvl_range, distance, nelements = None)
        data = SDECData.from_simulation(simulation_simple, packets_mode)
        
        # pytest.raises(ValueError, plotter._calculate_plotting_data, distance = -50 * u.Mpc, packets_mode = "virtual", packet_wvl_range = None, nelements = None)
        
        plot_frequency_bins = data.spectrum_frequency_bins # this changed
        plot_wavelength = data.spectrum_wavelength
        plot_freqency = data.spectrum_frequency
        
        # should I test SDECData here?
        if packet_wvl_range:
            packet_nu_range = packet_wvl_range.to("Hz", u.spectral())
            
            test_start_idx = (
                np.argmax(plot_frequency_bins > packet_nu_range[1]) - 1
            )
            test_end_idx = np.argmin(plot_frequency_bins < packet_nu_range[0])
            
            packet_wvl_range_mask = np.zeros(plot_wavelength.size, dtype=bool)
            packet_wvl_range_mask[test_start_idx:test_end_idx] = True

            assert (plotter.packet_wvl_range_mask[test_start_idx:test_end_idx]).all()
            np.testing.assert_allclose(plotter.plot_frequency_bins, plot_frequency_bins[test_start_idx: test_end_idx+1])
            np.testing.assert_allclose(plotter.packet_wvl_range_mask, packet_wvl_range_mask)
            np.testing.assert_allclose(plotter.plot_wavelength, plot_wavelength[packet_wvl_range_mask])
            np.testing.assert_allclose(plotter.plot_frequency, plot_freqency[packet_wvl_range_mask])
        else:
            np.testing.assert_allclose(plotter.packet_wvl_range_mask, np.ones(plot_wavelength.size, dtype = bool))
        
        if distance is None:
            assert plotter.lum_to_flux == 1
        else:
            if distance<=0:
                pytest.raises(ValueError)
            else:
                assert plotter.lum_to_flux == 4.0 * np.pi * (distance.to("cm")) ** 2
        
        emission_luminosities_df, emission_species = plotter._calculate_emission_luminosities(packets_mode=packets_mode, packet_wvl_range=packet_wvl_range)
        # pd.testing.assert_frame_equal(plotter.emission_luminosities_df, emission_luminosities_df)
        np.testing.assert_allclose(plotter.emission_species, emission_species)         
        
        absorption_luminosities_df, absorption_species = plotter._calculate_absorption_luminosities(packets_mode=packets_mode, packet_wvl_range=packet_wvl_range)
        np.testing.assert_allclose(plotter.absorption_species, absorption_species) 
        
        
        total_luminosities_df = (
            absorption_luminosities_df + emission_luminosities_df.drop(["noint", "escatter"], axis=1)
        )
        
        if nelements is None and plotter._species_list is None:
            np.testing.assert_allclose(plotter.species, np.array(list(total_luminosities_df.keys())))
        # elif plotter._species_list is not None:
            
        
        
        
            
            
            
            
            
            
            
            
        
        
        
    
    
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
                                
            if trace.name == "Other Elements":
                assert (trace.x == plotter.emission_luminosities_df.index.values).all()
                assert (trace.y == plotter.emission_luminosities_df.other.values).all()
        
        
    # TODO: Call generate_plot_mpl() with several PnCs of parameters and test
    # almost plotter's properties (esp. those saved by calculate_plotting_data)
    # against saved test data for those combinations (need to figure out a good
    # structure for saving test data for different PnCs)
