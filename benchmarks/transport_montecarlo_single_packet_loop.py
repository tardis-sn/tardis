"""
Basic TARDIS Benchmark.
"""

from numba.np.ufunc.parallel import get_num_threads, get_thread_id

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo import single_packet_loop


class BenchmarkTransportMontecarloSinglePacketLoop(BenchmarkBase):
    """
    Class to benchmark the single packet loop function.
    """

    def setup(self):
        self.Packet = self.packet
        self.geometry_1d = self.verysimple_numba_radial_1d_geometry
        self.time_explosion = self.verysimple_time_explosion
        self.opacity_state = self.verysimple_opacity_state
        self.estimator = (
            self.transport_state.radfield_mc_estimators.create_estimator_list(
                get_num_threads()
            )[get_thread_id()]
        )
        self.packet_collection_3d = self.verysimple_3vpacket_collection
        self.RpacketTracker = self.rpacket_tracker

    def time_single_packet_loop(self):
        single_packet_loop.single_packet_loop(
            self.Packet,
            self.geometry_1d,
            self.time_explosion,
            self.opacity_state,
            self.estimator,
            self.packet_collection_3d,
            self.RpacketTracker,
            self.montecarlo_configuration,
        )
