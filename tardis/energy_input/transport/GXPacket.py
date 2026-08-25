from enum import IntEnum

import numba as nb
import numpy as np
from numba import float64, int64
from numba.experimental import jitclass


class GXPacketStatus(IntEnum):
    """Status codes for gamma-ray packet transport."""

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
    ("time_idx", int64),
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
        time_idx,
        initialize_tau=True,
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
        self.time_idx = time_idx
        # TODO: rename to tau_event
        if initialize_tau:
            self.tau = -np.log(np.random.random())
        else:
            self.tau = 0.0

    def get_location_r(self):
        """Calculate radius of the packet

        Returns
        -------
            float: packet radius
        """
        return np.sqrt(
            self.location[0] ** 2.0 + self.location[1] ** 2.0 + self.location[2] ** 2.0
        )


@jitclass
class GXPacketCollection:
    """
    Gamma-ray packet collection
    """

    location: nb.float64[:, :]  # type: ignore[misc]
    direction: nb.float64[:, :]  # type: ignore[misc]
    energy_rf: nb.float64[:]  # type: ignore[misc]
    energy_cmf: nb.float64[:]  # type: ignore[misc]
    nu_rf: nb.float64[:]  # type: ignore[misc]
    nu_cmf: nb.float64[:]  # type: ignore[misc]
    status: nb.int64[:]  # type: ignore[misc]
    shell: nb.int64[:]  # type: ignore[misc]
    time_start: nb.float64[:]  # type: ignore[misc]
    time_index: nb.int64[:]  # type: ignore[misc]
    source_isotopes: nb.types.UnicodeCharSeq(16)[:]  # type: ignore[misc]
    packet_seeds: nb.int64[:]  # type: ignore[misc]

    def __init__(
        self,
        location: np.ndarray,
        direction: np.ndarray,
        energy_rf: np.ndarray,
        energy_cmf: np.ndarray,
        nu_rf: np.ndarray,
        nu_cmf: np.ndarray,
        status: np.ndarray,
        shell: np.ndarray,
        time_start: np.ndarray,
        time_index: np.ndarray,
        source_isotopes: np.ndarray,
        packet_seeds: np.ndarray,
    ) -> None:
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
        self.source_isotopes = source_isotopes
        self.packet_seeds = packet_seeds
