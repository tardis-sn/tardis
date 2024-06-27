"""
Basic TARDIS Benchmark.
"""

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo import single_packet_loop
from asv_runner.benchmarks.mark import skip_benchmark


class BenchmarkMontecarloMontecarloNumbaVpacket(BenchmarkBase):
    """
    Class to benchmark the single packet loop function.
    """

    @skip_benchmark
    def time_single_packet_loop(self):
        single_packet_loop.single_packet_loop(
            self.packet,
            self.verysimple_numba_radial_1d_geometry,
            self.verysimple_time_explosion,
            self.verysimple_opacity_state,
            self.estimators,
            self.verysimple_packet_collection,
            self.rpacket_tracker,
            self.montecarlo_configuration
        )
