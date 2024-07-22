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

    def time_run_tardis_rpacket_tracking(self):
        run_tardis(
            self.config_rpacket_tracking,
            atom_data=self.atomic_dataset,
            show_convergence_plots=False,
        )
