from numba import float64, int64, njit, from_dtype
from numba.experimental import jitclass
from numba.typed import List
import numpy as np
import pandas as pd


boundary_interaction_dtype = np.dtype(
    [
        ("event_id", "int64"),
        ("current_shell_id", "int64"),
        ("next_shell_id", "int64"),
    ]
)


line_interaction_dtype = np.dtype(
    [
        ("interaction_id", "int64"),
        ("shell_id", "int64"),
        ("r", "float64"),
        ("in_nu", "float64"),
        ("in_id", "int64"),
        ("out_id", "int64"),
    ]
)

rpacket_tracker_spec = [
    ("seed", int64),
    ("index", int64),
    ("status", int64[:]),
    ("r", float64[:]),
    ("nu", float64[:]),
    ("mu", float64[:]),
    ("energy", float64[:]),
    ("shell_id", int64[:]),
    ("interaction_type", int64[:]),
    ("boundary_interaction", from_dtype(boundary_interaction_dtype)[:]),
    ("line_interaction", from_dtype(line_interaction_dtype)[:]),
    ("num_interactions", int64),
    ("boundary_interactions_index", int64),
    ("line_interactions_index", int64),
    ("event_id", int64),
    ("extend_factor", int64),
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

    def __init__(self, length):
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
        self.interaction_type = np.empty(length, dtype=np.int64)
        self.boundary_interaction = np.empty(
            length,
            dtype=boundary_interaction_dtype,
        )
        self.line_interaction = np.empty(
            length,
            dtype=line_interaction_dtype,
        )
        self.num_interactions = 0
        self.boundary_interactions_index = 0
        self.line_interactions_index = 0
        self.event_id = 1
        self.extend_factor = 2

    def extend_array(self, array, array_length):
        temp_array = np.empty(
            array_length * self.extend_factor, dtype=array.dtype
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
            self.interaction_type = self.extend_array(
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

    def track_line_interaction(self, r_packet):
        """
        Track line interaction properties
        """
        if self.line_interactions_index >= self.line_interaction.size:
            self.line_interaction = self.extend_array(
                self.line_interaction,
                self.line_interaction.size,
            )

        self.line_interaction[self.line_interactions_index][
            "event_id"
        ] = self.event_id
        self.event_id += 1

        self.line_interaction[self.line_interactions_index][
            "shell_id"
        ] = r_packet.current_shell_id
        self.line_interaction[self.line_interactions_index]["r"] = r_packet.r
        self.line_interaction[self.line_interactions_index][
            "in_nu"
        ] = r_packet.last_interaction_in_nu
        self.line_interaction[self.line_interactions_index][
            "in_id"
        ] = r_packet.last_line_interaction_in_id
        self.line_interaction[self.line_interactions_index][
            "out_id"
        ] = r_packet.last_line_interaction_out_id
        self.line_interactions_index += 1


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
        self.line_interaction = self.line_interaction[
            : self.ine_interactions_index
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
    len_df = sum([len(tracker.r) for tracker in rpacket_trackers])
    index_array = np.empty([2, len_df], dtype="int")
    df_dtypes = np.dtype(
        [
            ("status", np.int64),
            ("seed", np.int64),
            ("r", np.float64),
            ("nu", np.float64),
            ("mu", np.float64),
            ("energy", np.float64),
            ("shell_id", np.int64),
            ("interaction_type", np.int64),
        ]
    )
    rpacket_tracker_ndarray = np.empty(len_df, df_dtypes)
    cur_index = 0
    for rpacket_tracker in rpacket_trackers:
        prev_index = cur_index
        cur_index = prev_index + len(rpacket_tracker.r)
        for j, column_name in enumerate(df_dtypes.fields.keys()):
            rpacket_tracker_ndarray[column_name][
                prev_index:cur_index
            ] = getattr(rpacket_tracker, column_name)
        index_array[0][prev_index:cur_index] = getattr(rpacket_tracker, "index")
        index_array[1][prev_index:cur_index] = range(cur_index - prev_index)
    return pd.DataFrame(
        rpacket_tracker_ndarray,
        index=pd.MultiIndex.from_arrays(index_array, names=["index", "step"]),
        columns=df_dtypes.names,
    )


rpacket_last_interaction_tracker_spec = [
    ("index", int64),
    ("r", float64),
    ("nu", float64),
    ("energy", float64),
    ("shell_id", int64),
    ("interaction_type", int64),
]


@jitclass(rpacket_last_interaction_tracker_spec)
class RPacketLastInteractionTracker(object):
    """
    Numba JITCLASS for storing the last interaction the RPacket undergoes.
    Parameters
    ----------
        index : int
            Index position of each RPacket
        r : float
            Radius of the shell where the RPacket is present
        nu : float
            Frequency of the RPacket
        energy : float
            Energy possessed by the RPacket
        shell_id : int
            Current Shell No in which the last interaction happened
        interaction_type: int
            Type of interaction the rpacket undergoes
    """

    def __init__(self):
        """
        Initialize properties with default values
        """
        self.index = -1
        self.r = -1.0
        self.nu = 0.0
        self.energy = 0.0
        self.shell_id = -1
        self.interaction_type = -1

    def track(self, r_packet):
        """
        Track properties of RPacket and override the previous values
        """
        self.index = r_packet.index
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = r_packet.last_interaction_type

    def finalize_array(self):
        """
        Added to make RPacketLastInteractionTracker compatible with RPacketTracker
        """
        pass

    def track_boundary_interaction(self, current_shell_id, next_shell_id):
        """
        Added to make RPacketLastInteractionTracker compatible with RPacketTracker
        """
        pass


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


@njit
def generate_rpacket_last_interaction_tracker_list(no_of_packets):
    """
    Parameters
    ----------
    no_of_packets : The count of RPackets that are sent in the ejecta

    Returns
    -------
    A list containing RPacketLastInteractionTracker for each RPacket
    """
    rpacket_trackers = List()
    for i in range(no_of_packets):
        rpacket_trackers.append(RPacketLastInteractionTracker())
    return rpacket_trackers
