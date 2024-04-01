import timeit
import os
from tardis.io.configuration.config_reader import Configuration
from tardis import run_tardis
from tardis.visualization import ConvergencePlots
from astropy import units as u

class Benchmarkvisualizationplot:
    """Class to benchmark TARDIS simulation and data analysis functions."""

    timeout = 200

    def setup(self):
        filename = "tardis_configv1_benchmark.yml"
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(dir_path, "data", filename)
        self.config = Configuration.from_yaml(path)
        self.config.atom_data = "kurucz_cd23_chianti_H_He.h5"
        self.sim = run_tardis(self.config, log_level="ERROR", show_progress_bars=False)

    def time_convergence_plot(self):
        plotter = ConvergencePlots.create_plasma_plot(self.sim)

    def benchmark_all(self):
        """Runs benchmarks for all functions and prints results."""
        print("Benchmarking TARDIS simulation...")
        setup_time = timeit.timeit(self.setup, number=1)
        simulation_time = timeit.timeit(self.time_run_tardis, number=1)
        print(f"Total setup time: {setup_time:.4f} seconds")
        print(f"Simulation time: {simulation_time:.4f} seconds")

        print("\nBenchmarking data analysis and visualization...")
        visualization_time = timeit.timeit(self.time_convergence_plot, number=1)
        print(f"Convergence plot time: {visualization_time:.4f} seconds")

if __name__ == "__main__":
    benchmark = Benchmarkvisualizationplot()
    benchmark.benchmark_all()
