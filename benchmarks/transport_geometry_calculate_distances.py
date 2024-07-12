from asv_runner.benchmarks.mark import parameterize

import tardis.transport.frame_transformations as frame_transformations
import tardis.transport.geometry.calculate_distances as calculate_distances
from benchmarks.benchmark_base import BenchmarkBase


class BenchmarkTransportGeometryCalculateDistances(BenchmarkBase):
    """
    Class to benchmark the calculate distances function.
    """

    @property
    def model(self):
        return 5.2e7

    def time_calculate_distance_boundary(self):
        mu = 0.3
        r = 7.5e14

        calculate_distances.calculate_distance_boundary(
            r, mu, self.geometry.r_inner[0], self.geometry.r_outer[0]
        )

    @parameterize(
        {
            "Parameters": [
                {
                    "packet": {
                        "nu_line": 0.1,
                        "is_last_line": True
                    },
                    "enable_full_relativity": True,
                },
                {
                    "packet": {
                        "nu_line": 0.2,
                        "is_last_line": False
                    },
                    "enable_full_relativity": True,
                }
            ]
        }
    )
    def time_calculate_distance_line(self, parameters):
        packet_params = parameters["packet"]
        nu_line = packet_params["nu_line"]
        is_last_line = packet_params["is_last_line"]
        enable_full_relativity = parameters["enable_full_relativity"]

        time_explosion = self.model

        doppler_factor = frame_transformations.get_doppler_factor(
            self.static_packet.r,
            self.static_packet.mu,
            time_explosion,
            enable_full_relativity
        )
        comov_nu = self.static_packet.nu * doppler_factor

        calculate_distances.calculate_distance_line(
            self.static_packet,
            comov_nu,
            is_last_line,
            nu_line,
            time_explosion,
            enable_full_relativity
        )
