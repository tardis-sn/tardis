"""
Basic TARDIS Benchmark.
"""
import functools

from asv_runner.benchmarks.mark import parameterize

from benchmarks.benchmark_base import BenchmarkBase
import tardis.transport.montecarlo.packets.trackers.tracker_full_util as tracker_full_util
import tardis.transport.montecarlo.packets.trackers.tracker_last_interaction_util as tracker_last_interaction_util


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
        tracker_full_util.trackers_full_to_df(self.TransportState.rpacket_tracker)

    def time_generate_rpacket_tracker_list(self):
        tracker_full_util.generate_tracker_full_list(50, 10)

    def time_generate_rpacket_last_interaction_tracker_list(self):
        tracker_last_interaction_util.generate_tracker_last_interaction_list(50)
