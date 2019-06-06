import numpy as np
from enum import IntEnum
from numba import int64, float64, boolean
from numba import jitclass, njit, gdb

from tardis.montecarlo.montecarlo_numba import njit_dict
from tardis import constants as const

class MonteCarloException(ValueError):
    pass

CLOSE_LINE_THRESHOLD = 1e-7
C_SPEED_OF_LIGHT = const.c.to('cm/s').value
MISS_DISTANCE = 1e99
SIGMA_THOMSON = const.sigma_T.to('cm^2').value
INVERSE_SIGMA_THOMSON = 1 / SIGMA_THOMSON

class PacketStatus(IntEnum):
    IN_PROCESS = 0
    EMITTED = 1
    REABSORBED = 2

class InteractionType(IntEnum):
    BOUNDARY = 1
    LINE = 2
    ESCATTERING = 3

rpacket_spec = [
    ('r', float64),
    ('mu', float64),
    ('nu', float64),
    ('energy', float64),
    ('next_line_id', int64),
    ('current_shell_id', int64),
    ('status', int64),
]

@njit(**njit_dict)
def calculate_distance_boundary(r, mu, r_inner, r_outer):
    delta_shell = 0
    if (mu > 0.0):
        # direction outward
        distance = np.sqrt(r_outer**2 + ((mu**2 - 1.0) * r**2)) - (r * mu)
        delta_shell = 1
    else:
        # going inward
        check = r_inner**2 + (r**2 * (mu**2 - 1.0))

        if (check >= 0.0):
            # hit inner boundary 
            distance = -r * mu - np.sqrt(check)
            delta_shell = -1
        else:
            # miss inner boundary 
            distance = np.sqrt(r_outer**2 + ((mu**2 - 1.0) * r**2)) - (r * mu)
            delta_shell = 1
    
    return distance, delta_shell

@njit(**njit_dict)
def calculate_distance_line(nu, comov_nu, nu_line, time_explosion):
    if nu_line == 0.0:
        return MISS_DISTANCE
    nu_diff = comov_nu - nu_line
    if np.abs(nu_diff / comov_nu) < CLOSE_LINE_THRESHOLD:
            nu_diff = 0.0
    if nu_diff >= 0:                    
        return (nu_diff / nu) * C_SPEED_OF_LIGHT * time_explosion
    else:
        print('nu difference is less than 0.0', nu_diff, comov_nu, nu, nu_line, time_explosion)
        raise MonteCarloException('nu difference is less than 0.0')

@njit(**njit_dict)
def calculate_distance_electron(electron_density, tau_event):
    return tau_event / (electron_density * SIGMA_THOMSON)

@njit(**njit_dict)
def calculate_tau_electron(electron_density, distance):    
    return electron_density * SIGMA_THOMSON * distance

@njit(**njit_dict)
def get_doppler_factor(r, mu, time_explosion):
    beta = (r / time_explosion) / C_SPEED_OF_LIGHT
    return 1.0 - mu * beta

@njit(**njit_dict)
def get_random_mu():
    return 2.0 * np.random.random() - 1.0




@jitclass(rpacket_spec)
class RPacket(object):
    def __init__(self, r, mu, nu, energy):
        self.r = r
        self.mu = mu
        self.nu = nu
        self.energy = energy
        self.current_shell_id = 0
        self.status = PacketStatus.IN_PROCESS

    def initialize_line_id(self, numba_plasma, numba_model):
        inverse_line_list_nu = numba_plasma.line_list_nu[::-1]
        doppler_factor = get_doppler_factor(self.r, self.mu,
                                            numba_model.time_explosion)
        comov_nu = self.nu * doppler_factor
        next_line_id = (len(numba_plasma.line_list_nu) -
                        np.searchsorted(inverse_line_list_nu, comov_nu))
        self.next_line_id = next_line_id



@njit(**njit_dict)
def trace_packet(r_packet, numba_model, numba_plasma):
    """

    Parameters
    ----------
    numba_model: tardis.montecarlo.montecarlo_numba.numba_interface.NumbaModel
    numba_plasma: tardis.montecarlo.montecarlo_numba.numba_interface.NumbaPlasma

    Returns
    -------

    """

    r_inner = numba_model.r_inner[r_packet.current_shell_id]
    r_outer = numba_model.r_outer[r_packet.current_shell_id]

    distance_boundary, delta_shell = calculate_distance_boundary(
        r_packet.r, r_packet.mu, r_inner, r_outer)

    # defining start for line interaction
    start_line_id = r_packet.next_line_id

    # defining taus
    tau_event = np.random.exponential()
    tau_trace_line_combined = 0.0

    # e scattering initialization

    cur_electron_density = numba_plasma.electron_density[
        r_packet.current_shell_id]
    distance_electron = calculate_distance_electron(
        cur_electron_density, tau_event)

    # Calculating doppler factor
    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu,
                                        numba_model.time_explosion)
    comov_nu = r_packet.nu * doppler_factor
    last_line = False
    cur_line_id = start_line_id
    for cur_line_id in range(start_line_id, len(numba_plasma.line_list_nu)):
        nu_line = numba_plasma.line_list_nu[cur_line_id]
        tau_trace_line = numba_plasma.tau_sobolev[
            cur_line_id, r_packet.current_shell_id]
            
        tau_trace_line_combined += tau_trace_line
        distance_trace = calculate_distance_line(
            r_packet.nu, comov_nu, nu_line, numba_model.time_explosion)
        tau_trace_electron = calculate_tau_electron(cur_electron_density,
                                                    distance_trace)
        tau_trace_combined = tau_trace_line_combined + tau_trace_electron

        if ((distance_boundary <= distance_trace) and
                (distance_boundary <= distance_electron)):
            interaction_type = InteractionType.BOUNDARY  # BOUNDARY
            r_packet.next_line_id = cur_line_id
            distance = distance_boundary
            break

        if ((distance_electron < distance_trace) and
                (distance_electron < distance_boundary)):
            interaction_type = InteractionType.ESCATTERING
            distance = distance_electron
            r_packet.next_line_id = cur_line_id
            break

        if tau_trace_combined > tau_event:
            interaction_type = InteractionType.LINE  # Line
            r_packet.next_line_id = cur_line_id
            distance = distance_trace
            break
    else:  # Executed when no break occurs in the for loop
        if cur_line_id == (len(numba_plasma.line_list_nu) - 1):
            # Treatment for last line
            cur_line_id += 1
        if distance_electron < distance_boundary:
            distance = distance_electron
            interaction_type = InteractionType.ESCATTERING
        else:
            distance = distance_boundary
            interaction_type = InteractionType.BOUNDARY

    r_packet.next_line_id = cur_line_id

    return distance, interaction_type, delta_shell


@njit(**njit_dict)
def move_rpacket(r_packet, distance, time_explosion, numba_estimator):
    """Move packet a distance and recalculate the new angle mu
    
    Parameters
    ----------
    distance : float
        distance in cm
    """


    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu, time_explosion)
    comov_nu = r_packet.nu * doppler_factor
    comov_energy = r_packet.energy * doppler_factor
    numba_estimator.j_estimator[r_packet.current_shell_id] += (
            comov_energy * distance)
    numba_estimator.nu_bar_estimator[r_packet.current_shell_id] += (
            comov_energy * distance * comov_nu)

    r = r_packet.r
    if (distance > 0.0):
        new_r = np.sqrt(r**2 + distance**2 +
                         2.0 * r * distance * r_packet.mu)
        r_packet.mu = (r_packet.mu * r + distance) / new_r
        r_packet.r = new_r

@njit(**njit_dict)
def move_packet_across_shell_boundary(packet, delta_shell,
                                      no_of_shells):
    """
    Move packet across shell boundary - realizing if we are still in the simulation or have
    moved out through the inner boundary or outer boundary and updating packet
    status.
    
    Parameters
    ----------
    distance : float
        distance to move to shell boundary
        
    delta_shell: int
        is +1 if moving outward or -1 if moving inward

    no_of_shells: int
        number of shells in TARDIS simulation
    """
    next_shell_id = packet.current_shell_id + delta_shell

    if next_shell_id >= no_of_shells:
        packet.status = PacketStatus.EMITTED
    elif next_shell_id < 0:
        packet.status = PacketStatus.REABSORBED
    else:
        packet.current_shell_id = next_shell_id

@njit(**njit_dict)
def scatter(r_packet, time_explosion):
    """
    General scattering for lines as well as thomson.
    2) get the doppler factor at that position with the old angle
    3) convert the current energy and nu into the comoving
        frame with the old mu
    4) Scatter and draw new mu - update mu
    5) Transform the comoving energy and nu back using the new mu
    
    Parameters
    ----------
    distance : [type]
        [description]
    """
    doppler_factor = get_doppler_factor(r_packet.r, r_packet.mu, time_explosion)
    comov_energy = r_packet.energy * doppler_factor
    comov_nu = r_packet.nu * doppler_factor
    r_packet.mu = get_random_mu()
    inverse_new_doppler_factor = 1. / get_doppler_factor(
        r_packet.r, r_packet.mu, time_explosion)
    r_packet.energy = comov_energy * inverse_new_doppler_factor
    r_packet.nu = comov_nu * inverse_new_doppler_factor
