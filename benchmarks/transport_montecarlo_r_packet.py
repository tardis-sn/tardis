"""
Basic TARDIS Benchmark.
"""
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.packet_trackers import (
    rpacket_trackers_to_dataframe,
)


class BenchmarkMontecarloMontecarloNumbaRPacket(BenchmarkBase):
    """
    Class to benchmark the numba R packet function.
    """

    def time_rpacket_trackers_to_dataframe(self):
        sim = self.simulation_rpacket_tracking_enabled
        transport_state = sim.transport.transport_state
        rpacket_trackers_to_dataframe(
            transport_state.rpacket_tracker
        )
