"""
Basic TARDIS Benchmark.
"""
import functools

from asv_runner.benchmarks.mark import parameterize

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo import packet_trackers


class BenchmarkTransportMontecarloPacketTrackers(BenchmarkBase):
    """
    Class to benchmark the numba R packet function.
    """

    repeat = 2

    @functools.cache
    def setup(self):
        sim = self.simulation_rpacket_tracking_enabled
        self.TransportState = sim.transport.transport_state

    def time_rpacket_trackers_to_dataframe(self):
        packet_trackers.rpacket_trackers_to_dataframe(
            self.TransportState.rpacket_tracker
        )

    def time_generate_rpacket_tracker_list(self):
        packet_trackers.generate_rpacket_tracker_list(10, 50)

    def time_generate_rpacket_last_interaction_tracker_list(self):
        packet_trackers.generate_rpacket_last_interaction_tracker_list(10)
