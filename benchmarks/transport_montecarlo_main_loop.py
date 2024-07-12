"""
Basic TARDIS Benchmark.
"""

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.montecarlo_main_loop import montecarlo_main_loop


class BenchmarkTransportMontecarloMainLoop(BenchmarkBase):
    """
    class to benchmark montecarlo_main_loop function.
    """

    def time_montecarlo_main_loop(self):
        montecarlo_main_loop(
            self.transport_state.packet_collection,
            self.transport_state.geometry_state,
            self.verysimple_time_explosion,
            self.transport_state.opacity_state,
            self.montecarlo_configuration,
            self.transport_state.radfield_mc_estimators,
            self.transport_state.spectrum_frequency.value,
            self.montecarlo_configuration.NUMBER_OF_VPACKETS,
            iteration=0,
            show_progress_bars=False,
            total_iterations=0
        )
