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

    repeat = 2

    def time_single_packet_loop(self):
        single_packet_loop.single_packet_loop(
            self.packet,
            self.verysimple_numba_radial_1d_geometry,
            self.verysimple_time_explosion,
            self.verysimple_opacity_state,
            self.transport_state.radfield_mc_estimators.create_estimator_list(
                get_num_threads()
            )[get_thread_id()],
            self.verysimple_3vpacket_collection,
            self.rpacket_tracker,
            self.montecarlo_configuration,
        )
