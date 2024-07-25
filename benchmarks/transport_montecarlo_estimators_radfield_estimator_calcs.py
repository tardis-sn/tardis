"""
Basic TARDIS Benchmark.
"""

import tardis.opacities.opacities as opacities
import tardis.transport.geometry.calculate_distances as calculate_distances
import tardis.transport.montecarlo.r_packet_transport as r_packet_transport
import tardis.transport.montecarlo.utils as utils
from benchmarks.benchmark_base import BenchmarkBase
from tardis.transport.montecarlo.estimators.radfield_estimator_calcs import (
    update_line_estimators,
)
from asv_runner.benchmarks.mark import parameterize


class BenchmarkMontecarloMontecarloNumbaPacket(BenchmarkBase):
    """
    Class to benchmark the numba packet function.
    """

    @parameterize(
        {
            "Parameters": [
                {
                    "cur_line_id": 0,
                    "distance_trace": 1e12,
                    "time_explosion": 5.2e7,
                    "enable_full_relativity": True,
                },
                {
                    "cur_line_id": 0,
                    "distance_trace": 0,
                    "time_explosion": 5.2e7,
                    "enable_full_relativity": True,
                },
                {
                    "cur_line_id": 1,
                    "distance_trace": 1e5,
                    "time_explosion": 1e10,
                    "enable_full_relativity": False,
                },
            ]
        }
    )
    def time_update_line_estimators(self, parameters):
        cur_line_id = parameters["cur_line_id"]
        distance_trace = parameters["distance_trace"]
        time_explosion = parameters["time_explosion"]
        enable_full_relativity = parameters["enable_full_relativity"]
        update_line_estimators(
            self.estimators,
            self.static_packet,
            cur_line_id,
            distance_trace,
            time_explosion,
            enable_full_relativity,
        )
