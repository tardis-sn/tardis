"""
Basic TARDIS Benchmark.
"""
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.packet_trackers import (
    rpacket_trackers_to_dataframe,
    generate_rpacket_tracker_list,
    generate_rpacket_last_interaction_tracker_list,
)


class BenchmarkTransportMontecarloPacketTrackers(BenchmarkBase):
    """
    Class to benchmark the numba R packet function.
    """

    def time_rpacket_trackers_to_dataframe(self):
        sim = self.simulation_rpacket_tracking_enabled
        transport_state = sim.transport.transport_state
        rpacket_trackers_to_dataframe(transport_state.rpacket_tracker)

    def time_generate_rpacket_tracker_list(self, no_of_packets, length):
        generate_rpacket_tracker_list(no_of_packets, length)

    def time_generate_rpacket_last_interaction_tracker_list(
        self, no_of_packets
    ):
        generate_rpacket_last_interaction_tracker_list(no_of_packets)

    time_generate_rpacket_tracker_list.params = ([1, 10, 50], [1, 10, 50])
    time_generate_rpacket_tracker_list.param_names = ["no_of_packets", "length"]

    time_generate_rpacket_last_interaction_tracker_list.params = [10, 100, 1000]
    time_generate_rpacket_last_interaction_tracker_list.param_names = [
        "no_of_packets"
    ]
