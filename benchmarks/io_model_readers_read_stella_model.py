"""
Basic TARDIS Benchmark.
"""
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.model import read_stella_model


# @skip_benchmark
class BenchmarkIoModelReadersReadStellaModel(BenchmarkBase):
    """
    Class to benchmark the run_tardis function.
    """

    def __init__(self):
        self.filename = "tardis/io/model/readers/tests/data/mesa.stella.dat"

    def time_read_stella_model(self):
        read_stella_model(self.get_absolute_path(self.filename))

    def time_read_stella_model_meta(self):
        stella_model_example_file = read_stella_model(self.get_absolute_path(self.filename))
        assert stella_model_example_file.metadata["zones"] == 400

    def time_read_stella_model_data(self):
        stella_model_example_file = read_stella_model(self.get_absolute_path(self.filename))
        assert stella_model_example_file.data.iloc[0, 0] == 6.006769337200000e29
        assert stella_model_example_file.data.iloc[-1, -1] == 2.123224906916000e04
