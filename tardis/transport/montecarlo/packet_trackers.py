from numba import float64, int64
from numba.experimental import jitclass
import numpy as np
import pandas as pd


#NOTE; in_nu, in_id, out_id corresponds to line interaction
#While using in_nu, in_id, out_id make sure that the corresponding...
#... interaction type is of line type only.
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
    ("in_nu", float64[:]),
    ("in_id", int64[:]),
    ("out_id", int64[:]),
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
        in_nu : float
            In case of line interaction, the frequency corresponding to the aborption process
        in_id : int
            In case of line interaction, the line id of the absorption process
        out_id : int
            In case of line interaction, the line id of the emission process
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
        self.in_nu = np.empty(self.length, dtype=np.float64)
        self.in_id = np.empty(self.length, dtype=np.int64)
        self.out_id = np.empty(self.length, dtype=np.int64)
        self.num_interactions = 0

    def track(self, r_packet_tracker):
        if self.num_interactions >= self.length:
            temp_length = self.length * 2
            temp_status = np.empty(temp_length, dtype=np.int64)
            temp_r = np.empty(temp_length, dtype=np.float64)
            temp_nu = np.empty(temp_length, dtype=np.float64)
            temp_mu = np.empty(temp_length, dtype=np.float64)
            temp_energy = np.empty(temp_length, dtype=np.float64)
            temp_shell_id = np.empty(temp_length, dtype=np.int64)
            temp_interaction_type = np.empty(temp_length, dtype=np.int64)
            temp_in_nu = np.empty(temp_length, dtype=np.float64)
            temp_in_id = np.empty(temp_length, dtype=np.int64)
            temp_out_id = np.empty(temp_length, dtype=np.int64)

            temp_status[: self.length] = self.status
            temp_r[: self.length] = self.r
            temp_nu[: self.length] = self.nu
            temp_mu[: self.length] = self.mu
            temp_energy[: self.length] = self.energy
            temp_shell_id[: self.length] = self.shell_id
            temp_interaction_type[: self.length] = self.interaction_type
            temp_in_nu[: self.length] = self.in_nu
            temp_in_id[: self.length] = self.in_id
            temp_out_id[: self.length] = self.out_id

            self.status = temp_status
            self.r = temp_r
            self.nu = temp_nu
            self.mu = temp_mu
            self.energy = temp_energy
            self.shell_id = temp_shell_id
            self.interaction_type = temp_interaction_type
            self.in_nu = temp_in_nu
            self.in_id = temp_in_id
            self.out_id = temp_out_id
            self.length = temp_length

        self.index = r_packet_tracker.index
        self.seed = r_packet_tracker.seed
        self.status[self.num_interactions] = r_packet_tracker.status
        self.r[self.num_interactions] = r_packet_tracker.r
        self.nu[self.num_interactions] = r_packet_tracker.nu
        self.mu[self.num_interactions] = r_packet_tracker.mu
        self.energy[self.num_interactions] = r_packet_tracker.energy
        self.shell_id[self.num_interactions] = r_packet_tracker.current_shell_id
        self.interaction_type[
            self.num_interactions
        ] = r_packet_tracker.last_interaction_type
        self.in_nu[self.num_interactions] = r_packet_tracker.last_interaction_in_nu
        self.in_id[self.num_interactions] = r_packet_tracker.last_line_interaction_in_id
        self.out_id[self.num_interactions] = r_packet_tracker.last_line_interaction_out_id
        self.num_interactions += 1

    def finalize_array(self):
        self.status = self.status[: self.num_interactions]
        self.r = self.r[: self.num_interactions]
        self.nu = self.nu[: self.num_interactions]
        self.mu = self.mu[: self.num_interactions]
        self.energy = self.energy[: self.num_interactions]
        self.shell_id = self.shell_id[: self.num_interactions]
        self.interaction_type = self.interaction_type[: self.num_interactions]
        self.in_nu = self.in_nu[: self.num_interactions]
        self.in_id = self.in_id[: self.num_interactions]
        self.out_id = self.out_id[: self.num_interactions]


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
    ("seed", int64),
    ("status", int64),
    ("r", float64),
    ("nu", float64),
    ("energy", float64),
    ("shell_id", int64),
    ("interaction_type", int64),
    ("in_nu", float64),
    ("in_id", int64),
    ("out_id", int64),
    ("line_r", float64),
]

#NOTE; in_nu, in_id, out_id, line_r doesn't necesarily correspond to the last interaction
#NOTE; they corresponds to the last line interaction
@jitclass(rpacket_last_interaction_tracker_spec)
class RPacketLastInteractionTracker(object):
    """
    Numba JITCLASS for storing the last interaction the RPacket undergoes.
    Parameters
    ----------
        index : int
            Index position of each RPacket
        seed : int
            Seed for each RPacket
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
        in_nu : float
            In case of line interaction, the frequency corresponding to the aborption process
        in_id : int
            In case of line interaction, the line id of the absorption process
        out_id : int
            In case of line interaction, the line id of the emission process
        line_r : float
            In case of line interaction, the radius at the last line interaction
    """

    def __init__(self):
        self.index = -1
        self.seed = 0
        self.status = 0
        self.r = -1.0
        self.nu = 0.0
        self.energy = 0.0
        self.shell_id = -1
        self.interaction_type = -1
        self.in_nu = 0.0
        self.in_id = -1
        self.out_id = -1
        self.line_r = 0.0

    def track(self, r_packet):
        self.index = r_packet.index
        self.seed = r_packet.seed
        self.status = r_packet.status
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = r_packet.last_interaction_type
        self.in_nu = r_packet.last_interaction_in_nu
        self.in_id = r_packet.last_line_interaction_in_id
        self.out_id = r_packet.last_line_interaction_out_id
        if self.interaction_type == 2:
            self.line_r = self.r
        
