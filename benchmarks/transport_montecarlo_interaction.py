"""
Basic TARDIS Benchmark.
"""

import tardis.transport.montecarlo.interaction as interaction
from benchmarks.benchmark_base import BenchmarkBase


# check it
class BenchmarkTransportMontecarloInteraction(BenchmarkBase):
    """
    Class to benchmark the numba interaction function.
    """

    def setup(self):
        self.Packet = self.packet
        self.time_explosion = self.verysimple_time_explosion
        self.enable_full_relativity = self.verysimple_enable_full_relativity
        self.opacity_state = self.verysimple_opacity_state

    def time_thomson_scatter(self):
        interaction.thomson_scatter(
            self.Packet, self.time_explosion, self.enable_full_relativity
        )

    def time_line_scatter(self):
        self.Packet.initialize_line_id(
            self.opacity_state, self.time_explosion, self.enable_full_relativity
        )
        interaction.line_scatter(
            self.Packet,
            self.time_explosion,
            "SCATTER",
            self.opacity_state,
            self.enable_full_relativity,
        )

    def time_line_emission(self):
        emission_line_id = 1000
        self.Packet.mu = 0.8599443103322428
        self.Packet.energy = 0.9114437898710559
        self.Packet.initialize_line_id(
            self.opacity_state, self.time_explosion, self.enable_full_relativity
        )

        interaction.line_emission(
            self.Packet,
            emission_line_id,
            self.time_explosion,
            self.opacity_state,
            self.enable_full_relativity,
        )
