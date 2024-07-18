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
                    "electron_density": 1e-5,
                    "tua_event": 1e10,
                },
                {
                    "electron_density": 1.0,
                    "tua_event": 1e10
                },
            ]
        }
    )
    def time_calculate_distance_electron(self, parameters):
        electron_density = parameters["electron_density"]
        tau_event = parameters["tua_event"]
        calculate_distances.calculate_distance_electron(
            electron_density, tau_event
        )

    @parameterize(
        {
            "Parameters": [
                {
                    "electron_density": 1e-5,
                    "distance": 1.0,
                },
                {
                    "electron_density": 1e10,
                    "distance": 1e10,
                },
                {
                    "electron_density": -1,
                    "distance": 0,
                },
                {
                    "electron_density": -1e10,
                    "distance": -1e10,
                },
            ]
        }
    )
    def time_calculate_tau_electron(self, parameters):
        electron_density = parameters["electron_density"]
        distance = parameters["distance"]
        opacities.calculate_tau_electron(electron_density, distance)

    def time_get_random_mu(self):
        utils.get_random_mu()

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

    @parameterize(
        {
            "Parameters": [
                {
                    "current_shell_id": 132,
                    "delta_shell": 11,
                    "no_of_shells": 132,
                },
                {
                    "current_shell_id": 132,
                    "delta_shell": 1,
                    "no_of_shells": 133,
                },
                {
                    "current_shell_id": 132,
                    "delta_shell": 2,
                    "no_of_shells": 133,
                },
            ]
        }
    )
    def time_move_packet_across_shell_boundary_emitted(self, parameters):
        current_shell_id = parameters["current_shell_id"]
        delta_shell = parameters["delta_shell"]
        no_of_shells = parameters["no_of_shells"]
        packet = self.packet
        packet.current_shell_id = current_shell_id
        r_packet_transport.move_packet_across_shell_boundary(
            packet, delta_shell, no_of_shells
        )

    @parameterize(
        {
            "Parameters": [
                {
                    "current_shell_id": 132,
                    "delta_shell": 132,
                    "no_of_shells": 132,
                },
                {
                    "current_shell_id": -133,
                    "delta_shell": -133,
                    "no_of_shells": -1e9,
                },
                {
                    "current_shell_id": 132,
                    "delta_shell": 133,
                    "no_of_shells": 133,
                },
            ]
        }
    )
    def time_move_packet_across_shell_boundary_reabsorbed(self, parameters):
        current_shell_id = parameters["current_shell_id"]
        delta_shell = parameters["delta_shell"]
        no_of_shells = parameters["no_of_shells"]
        packet = self.packet
        packet.current_shell_id = current_shell_id
        r_packet_transport.move_packet_across_shell_boundary(
            packet, delta_shell, no_of_shells
        )

    @parameterize(
        {
            "Parameters": [
                {
                    "current_shell_id": 132,
                    "delta_shell": -1,
                    "no_of_shells": 199,
                },
                {
                    "current_shell_id": 132,
                    "delta_shell": 0,
                    "no_of_shells": 132,
                },
                {
                    "current_shell_id": 132,
                    "delta_shell": 20,
                    "no_of_shells": 154,
                },
            ]
        }
    )
    def time_move_packet_across_shell_boundary_increment(self, parameters):
        current_shell_id = parameters["current_shell_id"]
        delta_shell = parameters["delta_shell"]
        no_of_shells = parameters["no_of_shells"]
        packet = self.packet
        packet.current_shell_id = current_shell_id
        r_packet_transport.move_packet_across_shell_boundary(
            packet, delta_shell, no_of_shells
        )
