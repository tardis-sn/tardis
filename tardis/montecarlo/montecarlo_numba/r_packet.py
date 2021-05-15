import math
from enum import IntEnum

import numpy as np
from numba import int64, float64, boolean
from numba import njit
from numba.experimental import jitclass

from tardis.montecarlo.montecarlo_numba import (
    njit_dict_no_parallel,
)

from tardis.montecarlo.montecarlo_numba.estimators import (
    set_estimators,
    set_estimators_full_relativity,
    update_line_estimators,
)
from tardis.montecarlo.montecarlo_numba.opacities import calculate_tau_electron
from tardis.montecarlo.montecarlo_numba.calculate_distances import (
    calculate_distance_boundary,
    calculate_distance_electron,
    calculate_distance_line,
)
from tardis.montecarlo.montecarlo_numba.frame_transformations import (
    get_doppler_factor,
)
from tardis.montecarlo import montecarlo_configuration
from tardis.montecarlo.montecarlo_numba.numba_config import (
    CLOSE_LINE_THRESHOLD,
    ENABLE_FULL_RELATIVITY,
)


class InteractionType(IntEnum):
    BOUNDARY = 1
    LINE = 2
    ESCATTERING = 3


class PacketStatus(IntEnum):
    IN_PROCESS = 0
    EMITTED = 1
    REABSORBED = 2


rpacket_spec = [
    ("r", float64),
    ("mu", float64),
    ("nu", float64),
    ("energy", float64),
    ("next_line_id", int64),
    ("current_shell_id", int64),
    ("status", int64),
    ("seed", int64),
    ("index", int64),
    ("is_close_line", boolean),
    ("last_interaction_type", int64),
    ("last_interaction_in_nu", float64),
    ("last_line_interaction_in_id", int64),
    ("last_line_interaction_out_id", int64),
]


@jitclass(rpacket_spec)
class RPacket(object):
    def __init__(self, r, mu, nu, energy, seed, index=0, is_close_line=False):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.current_shell_id = 0
        self.status = PacketStatus.IN_PROCESS
        self.seed = seed
        self.index = index
        self.is_close_line = is_close_line
        self.last_interaction_type = -1
        self.last_interaction_in_nu = 0.0
        self.last_line_interaction_in_id = -1
        self.last_line_interaction_out_id = -1

    def initialize_line_id(self, numba_plasma, numba_model):
        inverse_line_list_nu = numba_plasma.line_list_nu[::-1]
        doppler_factor = get_doppler_factor(
            self.r, self.mu, numba_model.time_explosion
        )
        comov_nu = self.nu * doppler_factor
        next_line_id = len(numba_plasma.line_list_nu) - np.searchsorted(
            inverse_line_list_nu, comov_nu
        )
        if next_line_id == len(numba_plasma.line_list_nu):
            next_line_id -= 1
        self.next_line_id = next_line_id


@njit(**njit_dict_no_parallel)
def trace_packet(r_packet, numba_model, numba_plasma, estimators):
    """
    Traces the RPacket through the ejecta and stops when an interaction happens (heart of the calculation)

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
    numba_model : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    numba_plasma : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma
    estimators : tardis.montecarlo.montecarlo_numba.numba_interface.Estimators

    Returns
    -------
    """

    r_inner = numba_model.r_inner[r_packet.current_shell_id]
    r_outer = numba_model.r_outer[r_packet.current_shell_id]

    (
        distance_boundary,
        delta_shell,
    ) = calculate_distance_boundary(r_packet.r, r_packet.mu, r_inner, r_outer)

    # defining start for line interaction
    start_line_id = r_packet.next_line_id

    # defining taus
    tau_event = -np.log(np.random.random())
    tau_trace_line_combined = 0.0

    # e scattering initialization

    cur_electron_density = numba_plasma.electron_density[
        r_packet.current_shell_id
    ]
    distance_electron = calculate_distance_electron(
        cur_electron_density, tau_event
    )

    # Calculating doppler factor
    doppler_factor = get_doppler_factor(
        r_packet.r, r_packet.mu, numba_model.time_explosion
    )
    comov_nu = r_packet.nu * doppler_factor

    cur_line_id = start_line_id  # initializing varibale for Numba
    # - do not remove
    last_line_id = len(numba_plasma.line_list_nu) - 1

    for cur_line_id in range(start_line_id, len(numba_plasma.line_list_nu)):

        # Going through the lines
        nu_line = numba_plasma.line_list_nu[cur_line_id]
        nu_line_last_interaction = numba_plasma.line_list_nu[cur_line_id - 1]

        # Getting the tau for the next line
        tau_trace_line = numba_plasma.tau_sobolev[
            cur_line_id, r_packet.current_shell_id
        ]

        # Adding it to the tau_trace_line_combined
        tau_trace_line_combined += tau_trace_line

        # Calculating the distance until the current photons co-moving nu
        # redshifts to the line frequency
        is_last_line = cur_line_id == last_line_id

        distance_trace = calculate_distance_line(
            r_packet,
            comov_nu,
            is_last_line,
            nu_line,
            numba_model.time_explosion,
        )

        # calculating the tau electron of how far the trace has progressed
        tau_trace_electron = calculate_tau_electron(
            cur_electron_density, distance_trace
        )

        # calculating the trace
        tau_trace_combined = tau_trace_line_combined + tau_trace_electron

        if (
            (distance_boundary <= distance_trace)
            and (distance_boundary <= distance_electron)
        ) and distance_trace != 0.0:
            interaction_type = InteractionType.BOUNDARY  # BOUNDARY
            r_packet.next_line_id = cur_line_id
            distance = distance_boundary
            break

        if (
            (distance_electron < distance_trace)
            and (distance_electron < distance_boundary)
        ) and distance_trace != 0.0:
            interaction_type = InteractionType.ESCATTERING
            # print('scattering')
            distance = distance_electron
            r_packet.next_line_id = cur_line_id
            break

        # Updating the J_b_lu and E_dot_lu
        # This means we are still looking for line interaction and have not
        # been kicked out of the path by boundary or electron interaction

        update_line_estimators(
            estimators,
            r_packet,
            cur_line_id,
            distance_trace,
            numba_model.time_explosion,
        )

        if (
            tau_trace_combined > tau_event
            and not montecarlo_configuration.disable_line_scattering
        ):
            interaction_type = InteractionType.LINE  # Line
            r_packet.last_interaction_in_nu = r_packet.nu
            r_packet.last_line_interaction_in_id = cur_line_id
            r_packet.next_line_id = cur_line_id
            distance = distance_trace
            break

        if not is_last_line:
            test_for_close_line(
                r_packet, cur_line_id + 1, nu_line, numba_plasma
            )

        # Recalculating distance_electron using tau_event -
        # tau_trace_line_combined
        distance_electron = calculate_distance_electron(
            cur_electron_density, tau_event - tau_trace_line_combined
        )

    else:  # Executed when no break occurs in the for loop
        # We are beyond the line list now and the only next thing is to see
        # if we are interacting with the boundary or electron scattering
        if cur_line_id == (len(numba_plasma.line_list_nu) - 1):
            # Treatment for last line
            cur_line_id += 1
        if distance_electron < distance_boundary:
            distance = distance_electron
            interaction_type = InteractionType.ESCATTERING
            # print('scattering')
        else:
            distance = distance_boundary
            interaction_type = InteractionType.BOUNDARY

    # r_packet.next_line_id = cur_line_id

    return distance, interaction_type, delta_shell


@njit(**njit_dict_no_parallel)
def move_r_packet(r_packet, distance, time_explosion, numba_estimator):
    """
    Move packet a distance and recalculate the new angle mu

    Parameters
    ----------
    r_packet : tardis.montecarlo.montecarlo_numba.r_packet.RPacket
        r_packet objects
    time_explosion : float
        time since explosion in s
    numba_estimator : tardis.montecarlo.montecarlo_numba.numba_interface.NumbaEstimator
        Estimators object
    distance : float
        distance in cm
    """

    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu, time_explosion)

    r = r_packet.r
    if distance > 0.0:
        new_r = np.sqrt(
            r * r + distance * distance + 2.0 * r * distance * r_packet.mu
        )
        r_packet.mu = (r_packet.mu * r + distance) / new_r
        r_packet.r = new_r

        comov_nu = r_packet.nu * doppler_factor
        comov_energy = r_packet.energy * doppler_factor

        if not ENABLE_FULL_RELATIVITY:
            set_estimators(
                r_packet, distance, numba_estimator, comov_nu, comov_energy
            )

        else:
            distance = distance * doppler_factor
            set_estimators_full_relativity(
                r_packet,
                distance,
                numba_estimator,
                comov_nu,
                comov_energy,
                doppler_factor,
            )


@njit(**njit_dict_no_parallel)
def move_packet_across_shell_boundary(packet, delta_shell, no_of_shells):
    """
    Move packet across shell boundary - realizing if we are still in the simulation or have
    moved out through the inner boundary or outer boundary and updating packet
    status.

    Parameters
    ----------
    distance : float
        distance to move to shell boundary

    delta_shell : int
        is +1 if moving outward or -1 if moving inward

    no_of_shells : int
        number of shells in TARDIS simulation
    """
    next_shell_id = packet.current_shell_id + delta_shell

    if next_shell_id >= no_of_shells:
        packet.status = PacketStatus.EMITTED
    elif next_shell_id < 0:
        packet.status = PacketStatus.REABSORBED
    else:
        packet.current_shell_id = next_shell_id


@njit(**njit_dict_no_parallel)
def test_for_close_line(r_packet, line_id, nu_line, numba_plasma):
    r_packet.is_close_line = abs(
        numba_plasma.line_list_nu[line_id] - nu_line
    ) < (nu_line * CLOSE_LINE_THRESHOLD)
