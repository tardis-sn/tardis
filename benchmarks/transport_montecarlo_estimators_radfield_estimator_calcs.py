"""
Basic TARDIS Benchmark.
"""

import functools

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_estimators_line,
)


class BenchmarkMontecarloMontecarloNumbaPacket(BenchmarkBase):
    """
    Class to benchmark the numba packet function.
    """

    repeat = 4

    @functools.cache
    def setup(self):
        self.EstimatorsLine = self.estimators_line
        self.StaticPacket = self.static_packet

    def time_update_estimators_line(self):
        update_estimators_line(
            self.EstimatorsLine, self.StaticPacket, 1, 1e12, 1e10, True
        )
