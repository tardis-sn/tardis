import numba as nb
import numpy as np
import pandas as pd
from numba import from_dtype, njit
from numba.experimental import jitclass
from numba.typed import List

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
        self.interaction_type = np.empty(length, dtype=np.int64)
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
        index_array[0][prev_index:cur_index] = rpacket_tracker.index
        index_array[1][prev_index:cur_index] = range(cur_index - prev_index)
    return pd.DataFrame(
        rpacket_tracker_ndarray,
        index=pd.MultiIndex.from_arrays(index_array, names=["index", "step"]),
        columns=df_dtypes.names,
    )


@jitclass
class RPacketLastInteractionTracker:
    index: nb.int64  # type: ignore[misc]
    r: nb.float64  # type: ignore[misc]
    nu: nb.float64  # type: ignore[misc]
    energy: nb.float64  # type: ignore[misc]
    shell_id: nb.int64  # type: ignore[misc]
    interaction_type: nb.int64  # type: ignore[misc]
    
    # LINE interaction tracking (before/after)
    line_before_nu: nb.float64  # type: ignore[misc]
    line_before_energy: nb.float64  # type: ignore[misc]
    line_before_mu: nb.float64  # type: ignore[misc]
    line_after_nu: nb.float64  # type: ignore[misc]
    line_after_energy: nb.float64  # type: ignore[misc]
    line_after_mu: nb.float64  # type: ignore[misc]

    # ESCATTERING interaction tracking (before/after)
    escat_before_nu: nb.float64  # type: ignore[misc]
    escat_before_energy: nb.float64  # type: ignore[misc]
    escat_before_mu: nb.float64  # type: ignore[misc]
    escat_after_nu: nb.float64  # type: ignore[misc]
    escat_after_energy: nb.float64  # type: ignore[misc]
    escat_after_mu: nb.float64  # type: ignore[misc]

    # CONTINUUM_PROCESS interaction tracking (before/after)
    continuum_before_nu: nb.float64  # type: ignore[misc]
    continuum_before_energy: nb.float64  # type: ignore[misc]
    continuum_before_mu: nb.float64  # type: ignore[misc]
    continuum_after_nu: nb.float64  # type: ignore[misc]
    continuum_after_energy: nb.float64  # type: ignore[misc]
    continuum_after_mu: nb.float64  # type: ignore[misc]    # BOUNDARY interaction tracking (single state)
    boundary_r: nb.float64  # type: ignore[misc]
    boundary_nu: nb.float64  # type: ignore[misc]
    boundary_energy: nb.float64  # type: ignore[misc]
    boundary_mu: nb.float64  # type: ignore[misc]
    boundary_from_shell_id: nb.int64  # type: ignore[misc]
    boundary_to_shell_id: nb.int64  # type: ignore[misc]
    
    # Counters for each interaction type
    line_interactions_count: nb.int64  # type: ignore[misc]
    escat_interactions_count: nb.int64  # type: ignore[misc]
    continuum_interactions_count: nb.int64  # type: ignore[misc]
    boundary_interactions_count: nb.int64  # type: ignore[misc]
    
    """
    Numba JITCLASS for storing the last interaction the RPacket undergoes,
    with specific tracking for each interaction type.

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
        line_before_* : float
            Properties before line interaction (nu, energy, mu)
        line_after_* : float
            Properties after line interaction (nu, energy, mu)
        escat_before_* : float
            Properties before electron scattering interaction (nu, energy, mu)
        escat_after_* : float
            Properties after electron scattering interaction (nu, energy, mu)
        continuum_before_* : float
            Properties before continuum process interaction (nu, energy, mu)
        continuum_after_* : float
            Properties after continuum process interaction (nu, energy, mu)
        boundary_* : float/int
            Properties during boundary crossing
        *_interactions_count : int
            Count of each interaction type
    """

    def __init__(self) -> None:
        """
        Initialize properties with default values
        """
        self.index = -1
        self.r = -1.0
        self.nu = 0.0
        self.energy = 0.0
        self.shell_id = -1
        self.interaction_type = -1

        # Initialize LINE interaction tracking
        self.line_before_nu = 0.0
        self.line_before_energy = 0.0
        self.line_before_mu = 0.0
        self.line_after_nu = 0.0
        self.line_after_energy = 0.0
        self.line_after_mu = 0.0

        # Initialize ESCATTERING interaction tracking
        self.escat_before_nu = 0.0
        self.escat_before_energy = 0.0
        self.escat_before_mu = 0.0
        self.escat_after_nu = 0.0
        self.escat_after_energy = 0.0
        self.escat_after_mu = 0.0

        # Initialize CONTINUUM_PROCESS interaction tracking
        self.continuum_before_nu = 0.0
        self.continuum_before_energy = 0.0
        self.continuum_before_mu = 0.0
        self.continuum_after_nu = 0.0
        self.continuum_after_energy = 0.0
        self.continuum_after_mu = 0.0

        # Initialize BOUNDARY interaction tracking
        self.boundary_r = -1.0
        self.boundary_nu = 0.0
        self.boundary_energy = 0.0
        self.boundary_mu = 0.0
        self.boundary_from_shell_id = -1
        self.boundary_to_shell_id = -1

        # Initialize counters
        self.line_interactions_count = 0
        self.escat_interactions_count = 0
        self.continuum_interactions_count = 0
        self.boundary_interactions_count = 0

    def track(self, r_packet) -> None:
        """
        Track properties of RPacket and override the previous values.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet to track.
        """
        self.index = r_packet.index
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = r_packet.last_interaction_type

    def track_line_interaction_before(self, r_packet) -> None:
        """
        Track packet state before line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before line interaction.
        """
        self.line_before_nu = r_packet.nu
        self.line_before_energy = r_packet.energy
        self.line_before_mu = r_packet.mu

    def track_line_interaction_after(self, r_packet) -> None:
        """
        Track packet state after line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after line interaction.
        """
        self.line_after_nu = r_packet.nu
        self.line_after_energy = r_packet.energy
        self.line_after_mu = r_packet.mu
        self.line_interactions_count += 1
        # Update general tracking
        self.track(r_packet)

    def track_escattering_interaction_before(self, r_packet) -> None:
        """
        Track packet state before electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before electron scattering.
        """
        self.escat_before_nu = r_packet.nu
        self.escat_before_energy = r_packet.energy
        self.escat_before_mu = r_packet.mu

    def track_escattering_interaction_after(self, r_packet) -> None:
        """
        Track packet state after electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after electron scattering.
        """
        self.escat_after_nu = r_packet.nu
        self.escat_after_energy = r_packet.energy
        self.escat_after_mu = r_packet.mu
        self.escat_interactions_count += 1
        # Update general tracking
        self.track(r_packet)

    def track_continuum_interaction_before(self, r_packet) -> None:
        """
        Track packet state before continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before continuum process.
        """
        self.continuum_before_nu = r_packet.nu
        self.continuum_before_energy = r_packet.energy
        self.continuum_before_mu = r_packet.mu

    def track_continuum_interaction_after(self, r_packet) -> None:
        """
        Track packet state after continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after continuum process.
        """
        self.continuum_after_nu = r_packet.nu
        self.continuum_after_energy = r_packet.energy
        self.continuum_after_mu = r_packet.mu
        self.continuum_interactions_count += 1
        # Update general tracking
        self.track(r_packet)

    def track_boundary_crossing(self, r_packet, from_shell_id, to_shell_id) -> None:
        """
        Track packet state during boundary crossing.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet during boundary crossing.
        from_shell_id : int
            Shell ID the packet is leaving.
        to_shell_id : int
            Shell ID the packet is entering.
        """
        self.boundary_r = r_packet.r
        self.boundary_nu = r_packet.nu
        self.boundary_energy = r_packet.energy
        self.boundary_mu = r_packet.mu
        self.boundary_from_shell_id = from_shell_id
        self.boundary_to_shell_id = to_shell_id
        self.boundary_interactions_count += 1
        # Update general tracking
        self.track(r_packet)

    def get_interaction_summary(self) -> nb.types.Tuple((nb.int64, nb.int64, nb.int64, nb.int64)):  # type: ignore[misc]
        """
        Get summary of all interaction counts.

        Returns
        -------
        tuple[int, int, int, int]
            Counts for (line, escattering, continuum, boundary) interactions.
        """
        return (
            self.line_interactions_count,
            self.escat_interactions_count,
            self.continuum_interactions_count,
            self.boundary_interactions_count,
        )

    def finalize_array(self):
        """
        Added to make RPacketLastInteractionTracker compatible with RPacketTracker
        """

    # Compatibility method for RPacketTracker interface
    def track_boundary_interaction(self, current_shell_id, next_shell_id):
        """
        Compatibility method for RPacketTracker interface.
        This method signature is maintained for backward compatibility,
        but actual boundary tracking should use the new method with r_packet.

        Parameters
        ----------
        current_shell_id : int
            Current shell ID
        next_shell_id : int
            Next shell ID
        """
        # Store minimal boundary information for compatibility
        self.boundary_from_shell_id = current_shell_id
        self.boundary_to_shell_id = next_shell_id
        self.boundary_interactions_count += 1


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
