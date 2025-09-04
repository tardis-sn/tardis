import pandas as pd

import numba as nb
import numpy as np
from numba import from_dtype, njit
from numba.experimental import jitclass
from numba.typed import List

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType

NO_INTERACTION_INT = int(InteractionType.NO_INTERACTION)

boundary_interaction_dtype = np.dtype(
    [
        ("event_id", "int64"),
        ("current_shell_id", "int64"),
        ("next_shell_id", "int64"),
    ]
)


@jitclass
class RPacketTracker:
    seed: nb.int64  # type: ignore[misc]
    index: nb.int64  # type: ignore[misc]
    status: nb.int64[:]  # type: ignore[misc]
    r: nb.float64[:]  # type: ignore[misc]
    nu: nb.float64[:]  # type: ignore[misc]
    mu: nb.float64[:]  # type: ignore[misc]
    energy: nb.float64[:]  # type: ignore[misc]
    shell_id: nb.int64[:]  # type: ignore[misc]
    interaction_type: nb.int64[:]  # type: ignore[misc]
    boundary_interaction: from_dtype(boundary_interaction_dtype)[:]  # type: ignore[misc]
    num_interactions: nb.int64  # type: ignore[misc]
    boundary_interactions_index: nb.int64  # type: ignore[misc]
    event_id: nb.int64  # type: ignore[misc]
    extend_factor: nb.int64  # type: ignore[misc]
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
            Frequency of the RPacket
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
        extend_factor : int
            The factor by which to extend the properties array when the size limit is reached
    """

    def __init__(self, length: int) -> None:
        """
        Initialize the variables with default value
        """
        self.seed = np.int64(0)
        self.index = np.int64(0)
        self.status = np.empty(length, dtype=np.int64)
        self.r = np.empty(length, dtype=np.float64)
        self.nu = np.empty(length, dtype=np.float64)
        self.mu = np.empty(length, dtype=np.float64)
        self.energy = np.empty(length, dtype=np.float64)
        self.shell_id = np.empty(length, dtype=np.int64)
        self.interaction_type = np.full(length, NO_INTERACTION_INT, dtype=np.int64)
        self.boundary_interaction = np.empty(
            length,
            dtype=boundary_interaction_dtype,
        )
        self.num_interactions = 0
        self.boundary_interactions_index = 0
        self.event_id = 1
        self.extend_factor = 2

    def extend_array(self, array, array_length):
        temp_array = np.empty(
            array_length * self.extend_factor, dtype=array.dtype
        )
        temp_array[:array_length] = array
        return temp_array

    def extend_interaction_type_array(self, array, array_length):
        temp_array = np.full(
            array_length * self.extend_factor, NO_INTERACTION_INT, dtype=array.dtype
        )
        temp_array[:array_length] = array
        return temp_array

    def track(self, r_packet):
        """
        Track important properties of RPacket
        """
        if self.num_interactions >= self.status.size:
            self.status = self.extend_array(self.status, self.status.size)
            self.r = self.extend_array(self.r, self.r.size)
            self.nu = self.extend_array(self.nu, self.nu.size)
            self.mu = self.extend_array(self.mu, self.mu.size)
            self.energy = self.extend_array(self.energy, self.energy.size)
            self.shell_id = self.extend_array(self.shell_id, self.shell_id.size)
            self.interaction_type = self.extend_interaction_type_array(
                self.interaction_type, self.interaction_type.size
            )

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

    def track_boundary_interaction(self, current_shell_id, next_shell_id):
        """
        Track boundary interaction properties
        """
        if self.boundary_interactions_index >= self.boundary_interaction.size:
            self.boundary_interaction = self.extend_array(
                self.boundary_interaction,
                self.boundary_interaction.size,
            )

        self.boundary_interaction[self.boundary_interactions_index][
            "event_id"
        ] = self.event_id
        self.event_id += 1

        self.boundary_interaction[self.boundary_interactions_index][
            "current_shell_id"
        ] = current_shell_id

        self.boundary_interaction[self.boundary_interactions_index][
            "next_shell_id"
        ] = next_shell_id

        self.boundary_interactions_index += 1

    def finalize_array(self):
        """
        Change the size of the array from length ( or multiple of length ) to
        the actual number of interactions
        """
        self.status = self.status[: self.num_interactions]
        self.r = self.r[: self.num_interactions]
        self.nu = self.nu[: self.num_interactions]
        self.mu = self.mu[: self.num_interactions]
        self.energy = self.energy[: self.num_interactions]
        self.shell_id = self.shell_id[: self.num_interactions]
        self.interaction_type = self.interaction_type[: self.num_interactions]
        self.boundary_interaction = self.boundary_interaction[
            : self.boundary_interactions_index
        ]


def rpacket_trackers_to_dataframe(rpacket_trackers):
    """Generates a dataframe from the rpacket_trackers list of RPacketCollection Objects.

    Parameters
    ----------
    rpacket_trackers : numba.typed.typedlist.List
        list of individual RPacketCollection class objects

    Returns
    -------
    pandas.core.frame.DataFrame
        Dataframe containing properties of RPackets as columns like status, seed, r, nu, mu, energy, shell_id, interaction_type

    """
    # Collect all data arrays by simply concatenating tracker arrays
    all_status = []
    all_seed = []
    all_r = []
    all_nu = []
    all_mu = []
    all_energy = []
    all_shell_id = []
    all_interaction_type = []
    all_packet_index = []
    all_step_index = []

    for i, rpacket_tracker in enumerate(rpacket_trackers):
        n_interactions = len(rpacket_tracker.r)

        all_status.append(rpacket_tracker.status)
        all_seed.append(np.full(n_interactions, rpacket_tracker.seed))
        all_r.append(rpacket_tracker.r)
        all_nu.append(rpacket_tracker.nu)
        all_mu.append(rpacket_tracker.mu)
        all_energy.append(rpacket_tracker.energy)
        all_shell_id.append(rpacket_tracker.shell_id)
        all_interaction_type.append(rpacket_tracker.interaction_type)
        all_packet_index.append(np.full(n_interactions, rpacket_tracker.index))
        all_step_index.append(np.arange(n_interactions))

    # Concatenate all arrays
    combined_interaction_type = np.concatenate(all_interaction_type)

    data = {
        "status": np.concatenate(all_status),
        "seed": np.concatenate(all_seed),
        "r": np.concatenate(all_r),
        "nu": np.concatenate(all_nu),
        "mu": np.concatenate(all_mu),
        "energy": np.concatenate(all_energy),
        "shell_id": np.concatenate(all_shell_id),
        "interaction_type": combined_interaction_type,
    }

    # Create multi-index
    index_arrays = [np.concatenate(all_packet_index), np.concatenate(all_step_index)]
    multi_index = pd.MultiIndex.from_arrays(index_arrays, names=["index", "step"])

    return pd.DataFrame(data, index=multi_index)


@njit
def generate_rpacket_tracker_list(no_of_packets, length):
    """
    Parameters
    ----------
    no_of_packets : The count of RPackets that are sent in the ejecta
    length : initial length of the tracking array

    Returns
    -------
    A list containing RPacketTracker for each RPacket
    """
    rpacket_trackers = List()
    for i in range(no_of_packets):
        rpacket_trackers.append(RPacketTracker(length))
    return rpacket_trackers