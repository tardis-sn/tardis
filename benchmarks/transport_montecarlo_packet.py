"""
Basic TARDIS Benchmark.
"""

import numpy as np
from asv_runner.benchmarks.mark import parameterize, skip_benchmark

import tardis.transport.montecarlo.estimators.radfield_mc_estimators
import tardis.transport.montecarlo.estimators.radfield_mc_estimators
import tardis.transport.montecarlo.numba_interface as numba_interface
import tardis.transport.montecarlo.opacities as opacities
import tardis.transport.montecarlo.r_packet as r_packet
import tardis.transport.montecarlo.utils as utils
import tardis.transport.frame_transformations as frame_transformations
import tardis.transport.geometry.calculate_distances as calculate_distances
import tardis.transport.r_packet_transport as r_packet_transport
from benchmarks.benchmark_base import BenchmarkBase
from tardis.model.geometry.radial1d import NumbaRadial1DGeometry
from tardis.transport.montecarlo.estimators.radfield_mc_estimators import (
    update_line_estimators,
)


class BenchmarkMontecarloMontecarloNumbaPacket(BenchmarkBase):
    """
    Class to benchmark the numba packet function.
    """

    @property
    def geometry(self):
        return NumbaRadial1DGeometry(
            r_inner=np.array([6.912e14, 8.64e14], dtype=np.float64),
            r_outer=np.array([8.64e14, 1.0368e15], dtype=np.float64),
            v_inner=np.array([-1, -1], dtype=np.float64),
            v_outer=np.array([-1, -1], dtype=np.float64),
        )

    @property
    def model(self):
        return numba_interface.NumbaModel(
            time_explosion=5.2e7,
        )

    @property
    def estimators(self):
        return tardis.transport.montecarlo.estimators.radfield_mc_estimators.RadiationFieldMCEstimators(
            j_estimator=np.array([0.0, 0.0], dtype=np.float64),
            nu_bar_estimator=np.array([0.0, 0.0], dtype=np.float64),
            j_blue_estimator=np.array(
                [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64
            ),
            Edotlu_estimator=np.array(
                [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], dtype=np.float64
            ),
            photo_ion_estimator=np.empty((0, 0), dtype=np.float64),
            stim_recomb_estimator=np.empty((0, 0), dtype=np.float64),
            bf_heating_estimator=np.empty((0, 0), dtype=np.float64),
            stim_recomb_cooling_estimator=np.empty((0, 0), dtype=np.float64),
            photo_ion_estimator_statistics=np.empty((0, 0), dtype=np.int64),
        )

    @parameterize(
        {
            "Packet params": [
                {"mu": 0.3, "r": 7.5e14},
                {"mu": -0.3, "r": 7.5e13},
                {"mu": -0.3, "r": 7.5e14},
            ]
        }
    )
    def time_calculate_distance_boundary(self, packet_params):
        mu = packet_params["mu"]
        r = packet_params["r"]

        calculate_distances.calculate_distance_boundary(
            r, mu, self.geometry.r_inner[0], self.geometry.r_outer[0]
        )

    @parameterize(
        {
            "Parameters": [
                {
                    "packet": {"nu_line": 0.1, "is_last_line": True},
                    "expected": None,
                },
                {
                    "packet": {"nu_line": 0.2, "is_last_line": False},
                    "expected": None,
                },
                {
                    "packet": {"nu_line": 0.5, "is_last_line": False},
                    "expected": utils.MonteCarloException,
                },
                {
                    "packet": {"nu_line": 0.6, "is_last_line": False},
                    "expected": utils.MonteCarloException,
                },
            ]
        }
    )
    def time_calculate_distance_line(self, parameters):
        packet_params = parameters["packet"]
        nu_line = packet_params["nu_line"]
        is_last_line = packet_params["is_last_line"]

        time_explosion = self.model.time_explosion

        doppler_factor = frame_transformations.get_doppler_factor(
            self.static_packet.r, self.static_packet.mu, time_explosion
        )
        comov_nu = self.static_packet.nu * doppler_factor

        calculate_distances.calculate_distance_line(
            self.static_packet,
            comov_nu,
            is_last_line,
            nu_line,
            time_explosion,
        )

    @parameterize(
        {
            "Parameters": [
                {
                    "electron_density": 1e-5,
                    "tua_event": 1e10,
                },
                {"electron_density": 1.0, "tua_event": 1e10},
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
                },
                {
                    "cur_line_id": 0,
                    "distance_trace": 0,
                    "time_explosion": 5.2e7,
                },
                {
                    "cur_line_id": 1,
                    "distance_trace": 1e5,
                    "time_explosion": 1e10,
                },
            ]
        }
    )
    def time_update_line_estimators(self, parameters):
        cur_line_id = parameters["cur_line_id"]
        distance_trace = parameters["distance_trace"]
        time_explosion = parameters["time_explosion"]
        update_line_estimators(
            self.estimators,
            self.static_packet,
            cur_line_id,
            distance_trace,
            time_explosion,
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

    @skip_benchmark
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
