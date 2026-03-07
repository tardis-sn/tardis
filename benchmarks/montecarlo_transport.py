"""
Basic TARDIS Benchmark.
"""

import functools

from asv_runner.benchmarks.mark import parameterize

from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.geometry.calculate_distances import (
    calculate_distance_boundary,
    calculate_distance_electron,
    calculate_distance_line,
)
from tardis.transport.montecarlo.packets.radiative_packet import RPacket


@parameterize({"Interaction type": ["scatter", "macroatom"]})
class BenchmarkMontecarloTransport(BenchmarkBase):
    """
    Class to benchmark the Monte Carlo transport module.
    """

    repeat = 2

    @functools.cache
    def setup(self, interaction_type):
        self.sim = self.nb_simulation_verysimple
        self.transport = self.sim.transport
        self.geometry = self.transport.transport_state.geometry_state
        self.opacity_state = self.sim.opacity_state
        self.macro_atom_state = self.sim.macro_atom_state
        self.time_explosion = self.transport.transport_state.time_explosion
        self.r_packet = RPacket(
            r=self.geometry.r_inner[0],
            mu=0.5,
            nu=1e15,
            energy=1.0,
            seed=1,
            index=0,
        )

    def time_initialize_transport_state(self, interaction_type):
        self.transport.initialize_transport_state(
            self.sim.simulation_state,
            self.opacity_state,
            self.macro_atom_state,
            self.sim.plasma,
            self.sim.no_of_packets,
            self.sim.no_of_virtual_packets,
            iteration=0,
        )

    def time_calculate_distance_boundary(self, interaction_type):
        calculate_distance_boundary(
            self.geometry.r_inner[0],
            0.5,
            self.geometry.r_inner[0],
            self.geometry.r_outer[0],
        )

    def time_calculate_distance_line(self, interaction_type):
        calculate_distance_line(
            self.r_packet,
            1e15,
            False,
            1e15,
            self.time_explosion,
            False,
        )

    def time_calculate_distance_electron(self, interaction_type):
        calculate_distance_electron(1e8, 1.0)
