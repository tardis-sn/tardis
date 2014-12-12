import numpy as np
cimport numpy as np

ctypedef np.int64_t int_type_t



cdef extern from "cmontecarlo.h":
    ctypedef enum rpacket_status_t:
        TARDIS_PACKET_STATUS_IN_PROCESS = 0
        TARDIS_PACKET_STATUS_EMITTED = 1
        TARDIS_PACKET_STATUS_REABSORBED = 2


    ctypedef struct rpacket_t:
        double nu
        double mu
        double energy
        double r
        double tau_event
        double nu_line
        int_type_t current_shell_id
        int_type_t next_line_id
        int_type_t last_line
        int_type_t close_line
        int_type_t recently_crossed_boundary
        int_type_t virtual_packet_flag
        int_type_t virtual_packet
        double d_line
        double d_electron
        double d_boundary
        rpacket_status_t next_shell_id


    double rpacket_get_nu (rpacket_t * packet)
    double rpacket_set_nu (rpacket_t * packet, double nu)


def create_rpacket():#
    cdef rpacket_t packet
    return packet
    # double nu, double mu, double energy)

cdef convert_dict2rpacket(packet_dict):
    cdef rpacket_t packet

    packet.nu = packet_dict['nu']
    packet.mu = packet_dict['mu']
    packet.energy = packet_dict['energy']
    packet.r = packet_dict['r']
    packet.tau_event = packet_dict['tau_event']
    packet.nu_line = packet_dict['nu_line']
    packet.current_shell_id = packet_dict['current_shell_id']
    packet.next_line_id = packet_dict['next_line_id']
    packet.last_line = packet_dict['last_line']
    packet.close_line = packet_dict['close_line']
    packet.recently_crossed_boundary = packet_dict['recently_crossed_boundary']
    packet.virtual_packet_flag = packet_dict['virtual_packet_flag']
    packet.virtual_packet = packet_dict['virtual_packet']
    packet.d_line = packet_dict['d_line']
    packet.d_electron = packet_dict['d_electron']
    packet.d_boundary = packet_dict['d_boundary']
    packet.next_shell_id = packet_dict['next_shell_id']

    return packet


def rpacket_get_nu_w(packet_dict):
    cdef rpacket_t packet
    packet = convert_dict2rpacket(packet_dict)
    return rpacket_get_nu(&packet)

def rpacket_set_nu_w(packet_dict, nu):
    cdef rpacket_t packet
    packet = convert_dict2rpacket(packet_dict)
    return rpacket_get_nu(&packet)

