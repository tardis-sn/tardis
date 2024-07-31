"""
Basic TARDIS Benchmark.
"""
from asv_runner.benchmarks.mark import parameterize

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo import packet_trackers


@parameterize({"num_packets": [10, 100], "length": [10, 50]})
class BenchmarkTransportMontecarloPacketTrackers(BenchmarkBase):
    """
    Class to benchmark the numba R packet function.
    """

    repeat = 2

    def setup(self, num_packets, length):
        sim = self.simulation_rpacket_tracking_enabled
        self.TransportState = sim.transport.transport_state

    def time_rpacket_trackers_to_dataframe(self, num_packets, length):
        packet_trackers.rpacket_trackers_to_dataframe(
            self.TransportState.rpacket_tracker
        )

    def time_generate_rpacket_tracker_list(self, num_packets, length):
        packet_trackers.generate_rpacket_tracker_list(num_packets, length)

    def time_generate_rpacket_last_interaction_tracker_list(
        self, num_packets, length
    ):
        packet_trackers.generate_rpacket_last_interaction_tracker_list(
            num_packets
        )
