"""
Basic TARDIS Benchmark.
"""
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.io.model.readers.generic_readers import read_simple_ascii_abundances, read_simple_ascii_density


# @skip_benchmark
class BenchmarkIoModelReadersAsciiReaders(BenchmarkBase):
    """
    Class to benchmark the ascii readers function.
    """

    def __init__(self):
        pass

    def time_simple_ascii_density_reader_time(self):
        read_simple_ascii_density(
            f"{self.example_model_file_dir}/tardis_simple_ascii_density_test.dat"
        )

    def time_simple_ascii_abundance_reader(self):
        read_simple_ascii_abundances(
            f"{self.example_model_file_dir}/artis_abundances.dat"
        )
