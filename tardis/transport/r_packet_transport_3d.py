import numpy as np
from numba import njit

from tardis.montecarlo.montecarlo_numba import njit_dict_no_parallel
from tardis.energy_input.util import (
    doppler_gamma
)


@njit(**njit_dict_no_parallel)
def move_packet_3d(packet, distance):
    """
    Moves packet a distance along its direction vector

    Parameters
    ----------
    packet : GXPacket object
    distance : float

    Returns
    -------
    packet : GXPacket object

    """
    location_old = packet.location
    direction = packet.direction

    location_new = location_old + distance * direction

    packet.location = location_new

    doppler_factor = doppler_gamma(
        packet.direction, packet.location, packet.time_current
    )

    packet.nu_cmf = packet.nu_rf * doppler_factor
    packet.energy_cmf = packet.energy_rf * doppler_factor

    return packet
