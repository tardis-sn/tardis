from tardis.base import run_tardis
import pytest
import pandas as pd
import numpy as np
from tardis.visualization.tools.sdec_plot import SDECData, SDECPlotter

class TestSDECPlotter:

    @pytest.fixture(scope="class")
    def plotter(self, config_verysimple, atomic_dataset):
        """Instantiate using a simple simulation"""
        # Setup simulation configuration using config_verysimple
        config_verysimple.montecarlo.iterations = 3
        config_verysimple.montecarlo.no_of_packets = 4000
        config_verysimple.montecarlo.last_no_of_packets = -1
        config_verysimple.spectrum.virtual.virtual_packet_logging = True

        atomic_data = deepcopy(atomic_dataset)
        sim = run_tardis(config_verysimple, atom_data=atomic_data)
        return SDECPlotter.from_simulation(sim)

    # Test each method of the instance's output with expected one (stored as df) or can test some aspects of df