"""
Basic TARDIS Benchmark.
"""

import numpy as np
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

import tardis.transport.montecarlo.interaction as interaction
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.numba_interface import (
    LineInteractionType,
)


@skip_benchmark
class BenchmarkMontecarloMontecarloNumbaInteraction(BenchmarkBase):
    """
    Class to benchmark the numba interaction function.
    """

    def time_thomson_scatter(self):
        packet = self.packet
        init_mu = packet.mu
        init_nu = packet.nu
        init_energy = packet.energy
        time_explosion = self.verysimple_time_explosion
        enable_full_relativity = self.verysimple_enable_full_relativity

        interaction.thomson_scatter(
            packet, time_explosion, enable_full_relativity
        )

        assert np.abs(packet.mu - init_mu) > 1e-7
        assert np.abs(packet.nu - init_nu) > 1e-7
        assert np.abs(packet.energy - init_energy) > 1e-7

    @parameterize(
        {
            "Line interaction type": [
                LineInteractionType.SCATTER,
                LineInteractionType.DOWNBRANCH,
                LineInteractionType.MACROATOM,
            ],
        }
    )
    def time_line_scatter(self, line_interaction_type):
        packet = self.packet
        init_mu = packet.mu
        init_nu = packet.nu
        init_energy = packet.energy
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
            self.verysimple_continuum_processes_enabled,
        )

        assert np.abs(packet.mu - init_mu) > 1e-7
        assert np.abs(packet.nu - init_nu) > 1e-7
        assert np.abs(packet.energy - init_energy) > 1e-7

    @parameterize(
        {
            "Test packet": [
                {
                    "mu": 0.8599443103322428,
                    "emission_line_id": 1000,
                    "energy": 0.9114437898710559,
                },
                {
                    "mu": -0.6975116557422458,
                    "emission_line_id": 2000,
                    "energy": 0.8803098648913266,
                },
                {
                    "mu": -0.7115661419975774,
                    "emission_line_id": 0,
                    "energy": 0.8800385929341252,
                },
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

        assert packet.next_line_id == emission_line_id + 1
