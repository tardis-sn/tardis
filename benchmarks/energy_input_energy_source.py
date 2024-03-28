"""
Basic TARDIS Benchmark.
"""
import numpy as np
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis.energy_input.samplers import create_energy_cdf


# @skip_benchmark
class BenchmarkEnergyInputEnergySource(BenchmarkBase):
    """
    Class to benchmark the run_tardis function.
    """

    @parameterize({"Energy CDF": [
            {
                "energy": np.array([100.0, 50.0]),
                "intensity": np.array([1.0, 1.0])
            },
            {
                "energy": np.array([50.0, 100.0]),
                "intensity": np.array([0.0, 1.0])
            },
    ]})
    def time_create_energy_cdf(self, values: dict):
        create_energy_cdf(values.get('energy'), values.get('intensity'))
