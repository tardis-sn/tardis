import numpy as np
cimport numpy as np

ctypedef np.int64_t int_type_t
ctypedef np.float64_t float_type_t



cdef extern from "cmontecarlo.h":
    ctypedef enum rpacket_status_t:
        TARDIS_PACKET_STATUS_IN_PROCESS = 0
        TARDIS_PACKET_STATUS_EMITTED = 1
        TARDIS_PACKET_STATUS_REABSORBED = 2


    ctypedef struct rpacket_t:
        pass

    double rpacket_get_nu (rpacket_t * packet)
    void rpacket_set_nu (rpacket_t * packet, double nu)

    double rpacket_get_mu (rpacket_t * packet)
    void rpacket_set_mu (rpacket_t * packet, double mu)

    double rpacket_get_energy (rpacket_t * packet)
    void rpacket_set_energy (rpacket_t * packet, double energy)

    double rpacket_get_r (rpacket_t * packet)
    void rpacket_set_r (rpacket_t * packet, double r)

    double rpacket_get_tau_event (rpacket_t * packet)
    void rpacket_set_tau_event (rpacket_t * packet, double tau_event)

    double rpacket_get_nu_line (rpacket_t * packet)
    void rpacket_set_nu_line (rpacket_t * packet, double nu_line)

    double rpacket_get_current_shell_id (rpacket_t * packet)
    void rpacket_set_current_shell_id (rpacket_t * packet, int_type_t current_shel_id)

    double rpacket_get_next_line_id (rpacket_t * packet)
    void rpacket_set_next_line_id (rpacket_t * packet, int_type_t next_line_id)

    double rpacket_get_last_line (rpacket_t * packet)
    void rpacket_set_last_line (rpacket_t * packet, int_type_t last_line)

    double rpacket_get_close_line (rpacket_t * packet)
    void rpacket_set_close_line (rpacket_t * packet, int_type_t close_line)

    double rpacket_get_recently_crossed_boundary (rpacket_t * packet)
    void rpacket_set_recently_crossed_boundary (rpacket_t * packet, int_type_t recently_crossed_boundary)

    double rpacket_get_virtual_packet_flag (rpacket_t * packet)
    void rpacket_set_virtual_packet_flag (rpacket_t * packet, int_type_t virtual_packet_flag)

    double rpacket_get_virtual_packet (rpacket_t * packet)
    void rpacket_set_virtual_packet (rpacket_t * packet, int_type_t virtual_packet)

    double rpacket_get_d_line (rpacket_t * packet)
    void rpacket_set_d_line (rpacket_t * packet, double d_line)

    double rpacket_get_d_electron (rpacket_t * packet)
    void rpacket_set_d_electron (rpacket_t * packet, double d_electron)

    double rpacket_get_d_boundary (rpacket_t * packet)
    void rpacket_set_d_boundary (rpacket_t * packet, double d_boundary)


cdef class RPacket:
    cdef rpacket_t * rpacket_ptr

    def __cinit__(self):
        cdef rpacket_t packet
        self.rpacket_ptr = &packet


    def __init__(self, rpacket_dict):
        for key in rpacket_dict:
            print key
            setattr(self, key, rpacket_dict[key])

    property nu:

        def __get__(self):
            return rpacket_get_nu(self.rpacket_ptr)

        def __set__(self, double nu):
            rpacket_set_nu(self.rpacket_ptr, nu)

    property mu:

        def __get__(self):
            return rpacket_get_mu(self.rpacket_ptr)

        def __set__(self, double mu):
            rpacket_set_mu(self.rpacket_ptr, mu)

    property energy:

        def __get__(self):
            return rpacket_get_energy(self.rpacket_ptr)

        def __set__(self, double energy):
            rpacket_set_energy(self.rpacket_ptr, energy)

    property r:

        def __get__(self):
            return rpacket_get_r(self.rpacket_ptr)

        def __set__(self, double r):
            rpacket_set_r(self.rpacket_ptr, r)

    property tau_event:

        def __get__(self):
            return rpacket_get_tau_event(self.rpacket_ptr)

        def __set__(self, double tau_event):
            rpacket_set_tau_event(self.rpacket_ptr, tau_event)

    property nu_line:
        def __get__(self):
            return rpacket_get_nu_line(self.rpacket_ptr)

        def __set__(self, double nu_line):
            rpacket_set_nu_line(self.rpacket_ptr, nu_line)

    property current_shell_id:
        def __get__(self):
            return rpacket_get_current_shell_id(self.rpacket_ptr)

        def __set__(self, int_type_t current_shell_id):
            rpacket_set_current_shell_id(self.rpacket_ptr, current_shell_id)

    property next_line_id:
        def __get__(self):
            return rpacket_get_next_line_id(self.rpacket_ptr)

        def __set__(self, int_type_t next_line_id):
            rpacket_set_next_line_id(self.rpacket_ptr, next_line_id)

    property last_line:
        def __get__(self):
            return rpacket_get_last_line(self.rpacket_ptr)

        def __set__(self, int_type_t last_line):
            rpacket_set_last_line(self.rpacket_ptr, last_line)

    property close_line:
        def __get__(self):
            return rpacket_get_close_line(self.rpacket_ptr)

        def __set__(self, int_type_t close_line):
            rpacket_set_close_line(self.rpacket_ptr, close_line)

    property recently_crossed_boundary:
        def __get__(self):
            return rpacket_get_recently_crossed_boundary(self.rpacket_ptr)

        def __set__(self, int_type_t recently_crossed_boundary):
            rpacket_set_recently_crossed_boundary(self.rpacket_ptr, recently_crossed_boundary)

    property virtual_packet_flag:
        def __get__(self):
            return rpacket_get_virtual_packet_flag(self.rpacket_ptr)

        def __set__(self, int_type_t virtual_packet_flag):
            rpacket_set_virtual_packet_flag(self.rpacket_ptr, virtual_packet_flag)

    property virtual_packet:
        def __get__(self):
            return rpacket_get_virtual_packet(self.rpacket_ptr)

        def __set__(self, int_type_t virtual_packet):
            rpacket_set_virtual_packet(self.rpacket_ptr, virtual_packet)

    property d_line:
        def __get__(self):
            return rpacket_get_d_line(self.rpacket_ptr)

        def __set__(self, double d_line):
            rpacket_set_d_line(self.rpacket_ptr, d_line)

    property d_electron:
        def __get__(self):
            return rpacket_get_d_electron(self.rpacket_ptr)

        def __set__(self, double d_electron):
            rpacket_set_d_electron(self.rpacket_ptr, d_electron)

    property d_boundary:
        def __get__(self):
            return rpacket_get_d_boundary(self.rpacket_ptr)

        def __set__(self, double d_boundary):
            rpacket_set_d_boundary(self.rpacket_ptr, d_boundary)

