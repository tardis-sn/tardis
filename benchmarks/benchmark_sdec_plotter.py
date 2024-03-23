import os
from tardis.io.configuration.config_reader import Configuration
from tardis import run_tardis
from tardis.visualization import SDECPlotter


class BenchmarkSDECplotter:
    """Class to benchmark Spectral element DEComposition (SDEC) Plot.
    """
    timeout = 200

    def setup(self):
        filename = "tardis_configv1_benchmark.yml"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(dir_path, "data", filename)
        config = Configuration.from_yaml(path)
        config.atom_data = "kurucz_cd23_chianti_H_He.h5"
        self.config = config
        self.sim = run_tardis(self.config, virtual_packet_logging=True)

    def time_sdec_plot_mpl(self):
        plotter = SDECPlotter.from_simulation(self.sim)
        plotter.generate_plot_mpl()

    def time_sdec_plot_ply(self):
        plotter = SDECPlotter.from_simulation(self.sim)
        plotter.generate_plot_ply()

