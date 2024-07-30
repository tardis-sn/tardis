"""
Basic TARDIS Benchmark.
"""
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo import packet_trackers


class BenchmarkTransportMontecarloPacketTrackers(BenchmarkBase):
    """
    Class to benchmark the numba R packet function.
    """

    def setup(self):
        sim = self.simulation_rpacket_tracking_enabled
        self.TransportState = sim.transport.transport_state

    def time_rpacket_trackers_to_dataframe(self):
        packet_trackers.rpacket_trackers_to_dataframe(
            self.TransportState.rpacket_tracker
        )

    def time_generate_rpacket_tracker_list(self, no_of_packets, length):
        packet_trackers.generate_rpacket_tracker_list(no_of_packets, length)

    time_generate_rpacket_tracker_list.params = ([1, 10, 50], [1, 10, 50])
    time_generate_rpacket_tracker_list.param_names = ["no_of_packets", "length"]

    def time_generate_rpacket_last_interaction_tracker_list(
        self, no_of_packets
    ):
        packet_trackers.generate_rpacket_last_interaction_tracker_list(
            no_of_packets
        )

    time_generate_rpacket_last_interaction_tracker_list.params = [10, 100, 1000]
    time_generate_rpacket_last_interaction_tracker_list.param_names = [
        "no_of_packets"
    ]
