import numpy as np
cimport numpy as np

ctypedef np.int64_t int_type_t



cdef extern from "cmontecarlo.h":
    ctypedef enum rpacket_status_t:
        TARDIS_PACKET_STATUS_IN_PROCESS = 0
        TARDIS_PACKET_STATUS_EMITTED = 1
        TARDIS_PACKET_STATUS_REABSORBED = 2


    ctypedef struct rpacket_t:
        pass

    double rpacket_get_nu (rpacket_t * packet)
    void rpacket_set_nu (rpacket_t * packet, double nu)


cdef class RPacket:
    cdef rpacket_t * rpacket_ptr

    def __cinit__(self):
        cdef rpacket_t packet
        self.rpacket_ptr = &packet

    property nu:

        def __get__(self):
            return rpacket_get_nu(self.rpacket_ptr)

        def __set__(self, double nu):
            print "HELLO"
            rpacket_set_nu(self.rpacket_ptr, nu)

