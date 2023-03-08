"""Benchmark the simulation."""
from tardis.simulation import Simulation
import os
from tardis.io.config_reader import Configuration

class BenchmarkSimulation:
    # params = ["tardis_configv1_verysimple.yml"]

    def setup(self):
        filename = "tardis_configv1_verysimple.yml"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(dir_path, filename)
        config = Configuration.from_yaml(path)
        config.atom_data = "/home/atharva/workspace/code/tardis-main/tardis/docs/kurucz_cd23_chianti_H_He.h5"
        sim = Simulation.from_config(config)
        self.config = config
        self.sim = sim

    def time_iterate(self):
        self.sim.iterate(int(self.config.montecarlo.no_of_packets))
    
    def time_advance_state(self):
        _ = self.sim.advance_state()

    
    




