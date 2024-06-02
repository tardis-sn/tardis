from numba import float64, int64
from numba.experimental import jitclass
import numpy as np


rpacket_tracker_spec = [
    ("length", int64),
    ("seed", int64),
    ("index", int64),
    ("status", int64[:]),
    ("r", float64[:]),
    ("nu", float64[:]),
    ("mu", float64[:]),
    ("energy", float64[:]),
    ("shell_id", int64[:]),
    ("interaction_type", int64[:]),
    ("num_interactions", int64),
]


@jitclass(rpacket_tracker_spec)
class RPacketTracker(object):
    """
    Numba JITCLASS for storing the information for each interaction a RPacket instance undergoes.
    Parameters
    ----------
        length : int
            Length of the initial array that is instantiated
        seed : int
            Seed for each RPacket
        index : int
            Index position of each RPacket
        status : int
            Current status of the RPacket as per interactions
        r : float
            Radius of the shell where the RPacket is present
        nu : float
            Luminosity of the RPacket
        mu : float
            Cosine of the angle made by the direction of movement of the RPacket from its original direction
        energy : float
            Energy possessed by the RPacket at a particular shell
        shell_id : int
            Current Shell No in which the RPacket is present
        interaction_type: int
            Type of interaction the rpacket undergoes
        num_interactions : int
            Internal counter for the interactions that a particular RPacket undergoes
    """

    def __init__(self, length):
        self.length = length
        self.seed = np.int64(0)
        self.index = np.int64(0)
        self.status = np.empty(self.length, dtype=np.int64)
        self.r = np.empty(self.length, dtype=np.float64)
        self.nu = np.empty(self.length, dtype=np.float64)
        self.mu = np.empty(self.length, dtype=np.float64)
        self.energy = np.empty(self.length, dtype=np.float64)
        self.shell_id = np.empty(self.length, dtype=np.int64)
        self.interaction_type = np.empty(self.length, dtype=np.int64)
        self.num_interactions = 0

    def track(self, r_packet):
        if self.num_interactions >= self.length:
            temp_length = self.length * 2
            temp_status = np.empty(temp_length, dtype=np.int64)
            temp_r = np.empty(temp_length, dtype=np.float64)
            temp_nu = np.empty(temp_length, dtype=np.float64)
            temp_mu = np.empty(temp_length, dtype=np.float64)
            temp_energy = np.empty(temp_length, dtype=np.float64)
            temp_shell_id = np.empty(temp_length, dtype=np.int64)
            temp_interaction_type = np.empty(temp_length, dtype=np.int64)

            temp_status[: self.length] = self.status
            temp_r[: self.length] = self.r
            temp_nu[: self.length] = self.nu
            temp_mu[: self.length] = self.mu
            temp_energy[: self.length] = self.energy
            temp_shell_id[: self.length] = self.shell_id
            temp_interaction_type[: self.length] = self.interaction_type

            self.status = temp_status
            self.r = temp_r
            self.nu = temp_nu
            self.mu = temp_mu
            self.energy = temp_energy
            self.shell_id = temp_shell_id
            self.interaction_type = temp_interaction_type
            self.length = temp_length

        self.index = r_packet.index
        self.seed = r_packet.seed
        self.status[self.num_interactions] = r_packet.status
        self.r[self.num_interactions] = r_packet.r
        self.nu[self.num_interactions] = r_packet.nu
        self.mu[self.num_interactions] = r_packet.mu
        self.energy[self.num_interactions] = r_packet.energy
        self.shell_id[self.num_interactions] = r_packet.current_shell_id
        self.interaction_type[
            self.num_interactions
        ] = r_packet.last_interaction_type
        self.num_interactions += 1

    def finalize_array(self):
        self.status = self.status[: self.num_interactions]
        self.r = self.r[: self.num_interactions]
        self.nu = self.nu[: self.num_interactions]
        self.mu = self.mu[: self.num_interactions]
        self.energy = self.energy[: self.num_interactions]
        self.shell_id = self.shell_id[: self.num_interactions]
        self.interaction_type = self.interaction_type[: self.num_interactions]


last_interaction_tracker_spec = [
    ("radius", float64[:]),
    ("shell_id", int64[:]),
    ("interaction_type", int64[:]),
    ("energy", float64[:]),
    ("in_id", int64[:]),
    ("in_nu", float64[:]),
    ("out_id", int64[:]),
    ("nu", float64[:]),
]


@jitclass(last_interaction_tracker_spec)
class RPacketLastInteractionTracker:
    def __init__(self, no_of_packets):
        self.radius = -1 * np.ones(no_of_packets, dtype=np.float64)
        self.shell_id = -1 * np.ones(no_of_packets, dtype=np.int32)
        self.nu = np.zeros(no_of_packets, dtype=np.float64)
        self.energy = np.zeros(no_of_packets, dtype=np.float64)
        self.interaction_type = -1 * np.ones(no_of_packets, dtype=np.int32)
        self.in_id = -1 * np.ones(no_of_packets, dtype=np.int32)
        self.in_nu = np.zeros(no_of_packets, dtype=np.float64)
        self.out_id = -1 * np.ones(no_of_packets, dtype=np.int32)

    def update_last_interaction(self, r_packet, i):
        self.radius[i] = r_packet.r
        self.shell_id[i] = r_packet.current_shell_id
        self.nu[i] = r_packet.nu
        self.energy[i] = r_packet.energy
        self.interaction_type[i] = r_packet.last_interaction_type
        self.in_id[i] = r_packet.last_line_interaction_in_id
        self.in_nu[i] = r_packet.last_interaction_in_nu
        self.out_id[i] = r_packet.last_line_interaction_out_id
