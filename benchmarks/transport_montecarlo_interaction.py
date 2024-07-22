"""
Basic TARDIS Benchmark.
"""

import tardis.transport.montecarlo.interaction as interaction
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.numba_interface import (
    LineInteractionType,
)
from asv_runner.benchmarks.mark import parameterize


class BenchmarkTransportMontecarloInteraction(BenchmarkBase):
    """
    Class to benchmark the numba interaction function.
    """

    def time_thomson_scatter(self):
        packet = self.packet
        time_explosion = self.verysimple_time_explosion
        enable_full_relativity = self.verysimple_enable_full_relativity

        interaction.thomson_scatter(
            packet, time_explosion, enable_full_relativity
        )

    @parameterize(
        {
            "Line interaction type": [
                LineInteractionType.SCATTER,
                LineInteractionType.MACROATOM,
            ],
        }
    )
    def time_line_scatter(self, line_interaction_type):
        packet = self.packet
        packet.initialize_line_id(
            self.verysimple_opacity_state,
            self.verysimple_time_explosion,
            self.verysimple_enable_full_relativity,
        )
        time_explosion = self.verysimple_time_explosion

        interaction.line_scatter(
            packet,
            time_explosion,
            line_interaction_type,
            self.verysimple_opacity_state,
            self.verysimple_enable_full_relativity,
        )

    @parameterize(
        {
            "Test packet": [
                {
                    "mu": 0.8599443103322428,
                    "emission_line_id": 1000,
                    "energy": 0.9114437898710559,
                }
            ]
        }
    )
    def time_line_emission(self, test_packet):
        emission_line_id = test_packet["emission_line_id"]
        packet = self.packet
        packet.mu = test_packet["mu"]
        packet.energy = test_packet["energy"]
        packet.initialize_line_id(
            self.verysimple_opacity_state,
            self.verysimple_time_explosion,
            self.verysimple_enable_full_relativity,
        )

        time_explosion = self.verysimple_time_explosion

        interaction.line_emission(
            packet,
            emission_line_id,
            time_explosion,
            self.verysimple_opacity_state,
            self.verysimple_enable_full_relativity,
        )
