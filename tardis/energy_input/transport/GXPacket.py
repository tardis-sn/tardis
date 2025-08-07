from enum import IntEnum

import numpy as np
from numba import float64, int64
from numba.experimental import jitclass

from tardis.energy_input.samplers import sample_decay_time, sample_energy
from tardis.energy_input.util import (
    H_CGS_KEV,
    doppler_factor_3d,
    get_index,
    get_random_unit_vector,
)


class GXPacketStatus(IntEnum):
    BETA_DECAY = -1
    COMPTON_SCATTER = 0
    PHOTOABSORPTION = 1
    PAIR_CREATION = 2
    IN_PROCESS = 3
    END = 4
    ESCAPED = 5


gxpacket_spec = [
    ("location", float64[:]),
    ("direction", float64[:]),
    ("energy_rf", float64),
    ("energy_cmf", float64),
    ("nu_rf", float64),
    ("nu_cmf", float64),
    ("status", int64),
    ("shell", int64),
    ("time_start", float64),
    ("time_index", int64),
    ("tau", float64),
]


@jitclass(gxpacket_spec)
class GXPacket:
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
        time_start,
        time_index,
    ):
        self.location = location
        self.direction = direction
        self.energy_rf = energy_rf
        self.energy_cmf = energy_cmf
        self.nu_rf = nu_rf
        self.nu_cmf = nu_cmf
        self.status = status
        self.shell = shell
        self.time_start = time_start
        self.time_index = time_index
        # TODO: rename to tau_event
        self.tau = -np.log(np.random.random())

    def get_location_r(self):
        """Calculate radius of the packet

        Returns
        -------
            float: packet radius
        """
        return np.sqrt(
            self.location[0] ** 2.0 + self.location[1] ** 2.0 + self.location[2] ** 2.0
        )


gxpacket_collection_spec = [
    ("location", float64[:, :]),
    ("direction", float64[:, :]),
    ("energy_rf", float64[:]),
    ("energy_cmf", float64[:]),
    ("nu_rf", float64[:]),
    ("nu_cmf", float64[:]),
    ("status", int64[:]),
    ("shell", int64[:]),
    ("time_start", float64[:]),
    ("time_index", int64[:]),
]

@jitclass(gxpacket_collection_spec)
class GXPacketCollection:
    """
    Gamma-ray packet collection
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
        time_start,
        time_index,
    ):
        self.location = location
        self.direction = direction
        self.energy_rf = energy_rf
        self.energy_cmf = energy_cmf
        self.nu_rf = nu_rf
        self.nu_cmf = nu_cmf
        self.status = status
        self.shell = shell
        self.time_start = time_start
        self.time_index = time_index
        #self.tau = -np.log(np.random.random())