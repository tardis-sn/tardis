"""
Basic TARDIS Benchmark.
"""
from benchmarks.benchmark_base import BenchmarkBase


class BenchmarkBasePlasma(BenchmarkBase):
    """
    Class to benchmark the BasePlasma solver factory in plasma/base.py.
    """

    repeat = 2

    def setup(self):
        """Set up the simulation and plasma for benchmarking."""
        self.sim = self.nb_simulation_verysimple
        self.plasma = self.sim.plasma
        self.j_blues = self.plasma.j_blues

    def time_build_graph(self):
        """Benchmark the time to build the plasma graph."""
        self.plasma._build_graph()

    def time_update_plasma(self):
        """Benchmark the time to update the plasma state via j_blues propagation."""
        self.plasma.update(j_blues=self.j_blues)

    def time_resolve_update_list(self):
        """Benchmark the time to resolve the update list for single variable."""
        self.plasma._resolve_update_list(["j_blues"])

    def time_get_value(self):
        """Benchmark the time to get the value of a single plasma property value."""
        self.plasma.get_value("t_rad")
