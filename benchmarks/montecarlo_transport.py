"""
Basic TARDIS Benchmark.
"""
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
    Class to benchmark core Monte Carlo transport distance calculations.
    """

    repeat = 2

    def setup(self, interaction_type):
        sim = self.nb_simulation_verysimple
        geometry = sim.transport.transport_state.geometry_state
        self.r_inner = geometry.r_inner[0]
        self.r_outer = geometry.r_outer[0]
        self.time_explosion = sim.transport.transport_state.time_explosion
        self.r_packet = RPacket(
            r=self.r_inner,
            mu=0.5,
            nu=1e15,
            energy=1.0,
            seed=1,
            index=0,
        )

    def time_calculate_distance_boundary(self, interaction_type):
        calculate_distance_boundary(
            self.r_inner,
            0.5,
            self.r_inner,
            self.r_outer,
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
