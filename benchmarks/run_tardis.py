"""
Basic TARDIS Benchmark.
"""

from benchmarks.benchmark_base import BenchmarkBase
from tardis import run_tardis

class BenchmarkRunTardis(BenchmarkBase):
    """
    Class to benchmark the `run tardis` function.
    """

    def __init__(self):
        super().__init__()
        filename = "data/tardis_configv1_benchmark.yml"
        self.path = self.get_relative_path(filename)

    def time_run_tardis(self):
        run_tardis(self.path, log_level="ERROR", show_progress_bars=False)
