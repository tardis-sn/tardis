from enum import IntEnum
import numpy as np
from numba import int64, float64
from numba.experimental import jitclass


class GXPacketStatus(IntEnum):
    BETA_DECAY = -1
    COMPTON_SCATTER = 0
    PHOTOABSORPTION = 1
    PAIR_CREATION = 2
    IN_PROCESS = 3
    END = 4


gxpacket_spec = [
    ("location", float64[:]),
    ("direction", float64[:]),
    ("energy_rf", float64),
    ("energy_cmf", float64),
    ("nu_rf", float64),
    ("nu_cmf", float64),
    ("status", int64),
    ("shell", int64),
    ("time_current", float64),
    ("tau", float64),
]


@jitclass(gxpacket_spec)
class GXPacket(object):
    """
    Indivisible gamma-ray packet
    """

    def __init__(
        self,
        location,
        direction,
        energy_rf,
        energy_cmf,
        nu_rf,
        nu_cmf,
        status,
        shell,
        time_current,
    ):
        self.location = location
        self.direction = direction
        self.energy_rf = energy_rf
        self.energy_cmf = energy_cmf
        self.nu_rf = nu_rf
        self.nu_cmf = nu_cmf
        self.status = status
        self.shell = shell
        self.time_current = time_current
        # TODO: rename to tau_event
        self.tau = -np.log(np.random.random())

    def get_location_r(self):
        """Calculate radius of the packet

        Returns:
            float: packet radius
        """
        return np.sqrt(
            self.location[0] ** 2.0
            + self.location[1] ** 2.0
            + self.location[2] ** 2.0
        )
