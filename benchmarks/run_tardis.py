"""
Basic TARDIS Benchmark.
"""

from benchmarks.benchmark_base import BenchmarkBase
from tardis import run_tardis

class BenchmarkRunTardis(BenchmarkBase):
    """
    Class to benchmark the `run tardis` function.
    """

    def time_run_tardis(self):
        run_tardis(
            self.config_verysimple,
            atom_data=self.atomic_dataset,
            show_convergence_plots=False,
        )
