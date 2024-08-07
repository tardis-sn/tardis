import functools

import tardis.transport.frame_transformations as frame_transformations
import tardis.transport.geometry.calculate_distances as calculate_distances
from benchmarks.benchmark_base import BenchmarkBase


class BenchmarkTransportGeometryCalculateDistances(BenchmarkBase):
    """
    Class to benchmark the calculate distances function.
    """

    @functools.cache
    def setup(self):
        self.StaticPacket = self.static_packet
        self.Geometry = self.geometry
        doppler_factor = frame_transformations.get_doppler_factor(
            self.StaticPacket.r, self.StaticPacket.mu, self.model, True
        )
        self.comov_nu = self.StaticPacket.nu * doppler_factor

    @property
    def model(self):
        return 5.2e7

    def time_calculate_distance_boundary(self):
        mu = 0.3
        r = 7.5e14

        calculate_distances.calculate_distance_boundary(
            r, mu, self.Geometry.r_inner[0], self.Geometry.r_outer[0]
        )

    def time_calculate_distance_line(self):
        nu_line = 0.2
        calculate_distances.calculate_distance_line(
            self.StaticPacket,
            self.comov_nu,
            True,
            nu_line,
            True,
            True,
        )
