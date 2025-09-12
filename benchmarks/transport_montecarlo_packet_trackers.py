"""
Basic TARDIS Benchmark.
"""
import functools

from asv_runner.benchmarks.mark import parameterize

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.packets import packet_trackers
import tardis.transport.montecarlo.packets.trackers.tracker_full


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
        tardis.transport.montecarlo.packets.trackers.tracker_full.trackers_full_to_df(
            self.TransportState.rpacket_tracker
        )

    def time_generate_rpacket_tracker_list(self):
        tardis.transport.montecarlo.packets.trackers.tracker_full.generate_tracker_full_list(50, 10)

    def time_generate_rpacket_last_interaction_tracker_list(self):
        packet_trackers.generate_rpacket_last_interaction_tracker_list(50)
