"""
Basic TARDIS Benchmark.
"""

from benchmarks.benchmark_base import BenchmarkBase
from tardis import run_tardis
from tardis.io.atom_data import AtomData


class BenchmarkRunTardis(BenchmarkBase):
    """
    Class to benchmark the `run tardis` function.
    """

    repeat = 2

    def setup(self):
        self.config = self.config_verysimple
        self.atom_data = AtomData.from_hdf(self.atomic_data_fname)

    def time_run_tardis(self):
        run_tardis(
            self.config, atom_data=self.atom_data, show_convergence_plots=False
        )
