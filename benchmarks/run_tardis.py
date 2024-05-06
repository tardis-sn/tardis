"""
Basic TARDIS Benchmark.
"""

from benchmarks.benchmark_base import BenchmarkBase
from tardis import run_tardis
from tardis.io.configuration.config_reader import Configuration


class BenchmarkRunTardis(BenchmarkBase):
    """
    Class to benchmark the `run tardis` function.
    """

    def __init__(self):
        super().__init__()
        self.config = None

    def setup(self):
        filename = "data/tardis_configv1_benchmark.yml"
        path = self.get_relative_path(filename)
        self.config = Configuration.from_yaml(path)

    def time_run_tardis(self):
        run_tardis(self.config, log_level="ERROR", show_progress_bars=False)
