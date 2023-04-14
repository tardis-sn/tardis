"""Basic TARDIS Benchmark."""
import os
from tardis.io.config_reader import Configuration
from tardis import run_tardis

class Benchmarkruntardis:
    """Class to benchmark the run_tardis function.
    """
    timeout = 200
    
    def setup(self):
        filename = "tardis_configv1_verysimple.yml"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(dir_path, "data", filename)
        config = Configuration.from_yaml(path)
        config.atom_data = "kurucz_cd23_chianti_H_He.h5"
        self.config = config
    
    def time_run_tardis(self):
        sim = run_tardis(self.config, log_level="ERROR", show_progress_bars=False)
