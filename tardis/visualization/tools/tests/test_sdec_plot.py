from tardis.base import run_tardis
import pytest
import pandas as pd
import numpy as np
from tardis.visualization.tools.sdec_plot import SDECData, SDECPlotter


class TestSDECPlotter:
    @pytest.fixture(scope="class")
    def plotter(self, config_verysimple, atomic_dataset):
        """Instantiate SDEC plotter using a simple simulation model"""
        # Setup simulation configuration using config_verysimple and
        # override properties in such a way to make the simulation run faster
        config_verysimple.montecarlo.iterations = 3
        config_verysimple.montecarlo.no_of_packets = 4000
        config_verysimple.montecarlo.last_no_of_packets = -1
        config_verysimple.spectrum.virtual.virtual_packet_logging = True
        config_verysimple.spectrum.num = 2000

        atomic_data = deepcopy(atomic_dataset)
        sim = run_tardis(config_verysimple, atom_data=atomic_data)
        return SDECPlotter.from_simulation(sim)

    # TODO: Call generate_plot_mpl() with several PnCs of parameters and test
    # almost plotter's properties (esp. those saved by calculate_plotting_data)
    # against saved test data for those combinations (need to figure out a good
    # structure for saving test data for different PnCs)
