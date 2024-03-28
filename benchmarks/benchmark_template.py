"""
Basic TARDIS Benchmark.
"""
from asv_runner.benchmarks.mark import skip_benchmark
from benchmarks.benchmark_base import BenchmarkBase


@skip_benchmark
class BenchmarkXx(BenchmarkBase):
    """
    Class to benchmark the Xx function.
    """

    def __init__(self):
        pass

    def time_template(self):
        pass
