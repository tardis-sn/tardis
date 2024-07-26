"""
Basic TARDIS Benchmark.
"""

import functools

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_line_estimators,
)


class BenchmarkMontecarloMontecarloNumbaPacket(BenchmarkBase):
    """
    Class to benchmark the numba packet function.
    """

    @functools.cache
    def setup(self):
        self.Estimators = self.estimators
        self.StaticPacket = self.static_packet

    def time_update_line_estimators(self):
        update_line_estimators(
            self.Estimators, self.StaticPacket, 1, 1e12, 1e10, True
        )
