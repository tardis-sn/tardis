"""
Basic TARDIS Benchmark.
"""
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.packet_trackers import (
    rpacket_trackers_to_dataframe,
)


class BenchmarkTransportMontecarloPacketTrackers(BenchmarkBase):

    """
    Class to benchmark the numba R packet function.
    """

    def setup(self):
        sim = self.simulation_rpacket_tracking_enabled
        self.TransportState = sim.transport.transport_state

    def time_rpacket_trackers_to_dataframe(self):
        rpacket_trackers_to_dataframe(self.TransportState.rpacket_tracker)
