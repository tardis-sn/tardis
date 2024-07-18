from enum import IntEnum
import numpy as np
from numba import int64, float64
from numba.experimental import jitclass

from tardis.energy_input.samplers import sample_decay_time, sample_energy
from tardis.energy_input.util import (
    get_index,
    get_random_unit_vector,
    doppler_factor_3d,
    H_CGS_KEV,
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
        self.number_of_packets = len(self.energy_rf)
        self.tau = -np.log(np.random.random(self.number_of_packets))


# @njit(**njit_dict_no_parallel)
def initialize_packet_properties(
    isotope_energy,
    isotope_intensity,
    positronium_energy,
    positronium_intensity,
    positronium_fraction,
    packet_energy,
    k,
    tau_start,
    tau_end,
    initial_radius,
    times,
    effective_times,
    average_power_per_mass,
    uniform_packet_energies=True,
):
    """Initialize the properties of an individual packet

    Parameters
    ----------
    isotope_energy : numpy array
        _description_
    isotope_intensity : numpy array
        _description_
    positronium_energy : numpy array
        _description_
    positronium_intensity : numpy array
        _description_
    positronium_fraction : float
        _description_
    packet_energy : float
        _description_
    k : int
        _description_
    tau_start : float
        _description_
    tau_end : float
        _description_
    initial_radius : float
        _description_
    times : numpy array
        _description_
    effective_times : numpy array
        _description_

    Returns
    -------
    _type_
        _description_
    """
    decay_time = np.inf

    # packets sampled uniformly in time
    if not uniform_packet_energies:
        z = np.random.random()
        decay_time = z * times[0] + (1 - z) * times[-1]
    else:
        decay_time = sample_decay_time(
            tau_start,
            end_tau=tau_end,
            decay_time_min=0,
            decay_time_max=times[-1],
        )

    decay_time_index = get_index(decay_time, times)

    cmf_energy = sample_energy(isotope_energy, isotope_intensity)

    z = np.random.random()
    if z < positronium_fraction:
        z = np.random.random()
        if cmf_energy == 511 and z > 0.25:
            cmf_energy = sample_energy(
                positronium_energy, positronium_intensity
            )

    energy_factor = 1.0
    if decay_time < times[0]:
        energy_factor = decay_time / times[0]
        decay_time = times[0]

    # compute the new energy factor based on the ratio of the end-chain decay power
    # to the average decay power
    if not uniform_packet_energies:
        energy_factor = (
            get_chain_decay_power_per_ejectamass(
                inventory,
                decay_time,
                times[0],
                isotope_energy * isotope_intensity / 100,
                positronium_energy * positronium_intensity / 100,
                tau_end,
            )
            / average_power_per_mass
        )

    scaled_r = initial_radius * effective_times[decay_time_index]

    initial_location = scaled_r * get_random_unit_vector()
    initial_direction = get_random_unit_vector()
    initial_energy = packet_energy * energy_factor

    # draw a random gamma-ray in shell
    packet = GXPacket(
        initial_location,
        initial_direction,
        1.0,
        initial_energy,
        0.0,
        0.0,
        GXPacketStatus.IN_PROCESS,
        k,
        decay_time,
    )

    packet.energy_rf = packet.energy_cmf / doppler_factor_3d(
        packet.direction,
        packet.location,
        packet.time_current,
    )

    packet.nu_cmf = cmf_energy / H_CGS_KEV

    packet.nu_rf = packet.nu_cmf / doppler_factor_3d(
        packet.direction,
        packet.location,
        packet.time_current,
    )

    return packet, decay_time_index
