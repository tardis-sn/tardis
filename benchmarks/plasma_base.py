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
        self.sim = self.nb_simulation_verysimple
        self.plasma = self.sim.plasma
        self.j_blues = self.plasma.j_blues

    def time_build_graph(self):
        self.plasma._build_graph()

    def time_update_plasma(self):
        self.plasma.update(j_blues=self.j_blues)

    def time_resolve_update_list(self):
        self.plasma._resolve_update_list(["j_blues"])

    def time_get_value(self):
        self.plasma.get_value("t_rad")
