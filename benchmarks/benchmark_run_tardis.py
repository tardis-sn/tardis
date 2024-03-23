"""Basic TARDIS Benchmark."""
import os
from tardis.io.configuration.config_reader import Configuration
from tardis import run_tardis
from astropy import units as u

class Benchmarkruntardis:
    """Class to benchmark the run_tardis function.
    """
    timeout = 200
    params = ['branch85_w7', 'uniform','power_law','exponential']
    param_names = ['Density']

    def setup(self, d):
        filename = "tardis_configv1_benchmark.yml"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(dir_path, "data", filename)
        config = Configuration.from_yaml(path)
        config.atom_data = "kurucz_cd23_chianti_H_He.h5"
        config.model.structure.density.type = d
        config.model.structure.velocity.start = 1000 * u.km/u.s
        config.model.structure.velocity.stop = 2000 * u.km/u.s
        config.model.structure.velocity.num = 20
        if d == 'uniform':
            config.model.structure.density.time_0 = 1 * u.day
            config.model.structure.density.value = 5e-10 * u.kg/u.cm**3
        elif d == 'power_law':
            config.model.structure.density.time_0 = 1 * u.day
            config.model.structure.density.rho_0 = 5e-10 * u.kg/u.cm**3
            config.model.structure.density.v_0 = 500 * u.km/u.s
            config.model.structure.density.exponent = -2
        elif d == 'exponential':
            config.model.structure.density.time_0 = 1 * u.day
            config.model.structure.density.rho_0 = 5e-10 * u.kg/u.cm**3
            config.model.structure.density.v_0 = 500 * u.km/u.s
        self.config = config

    def time_run_tardis(self,d):
        sim = run_tardis(self.config, log_level="ERROR",
                         show_progress_bars=False)
