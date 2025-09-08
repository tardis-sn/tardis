import numba as nb
import numpy as np
import pandas as pd
from numba import njit
from numba.experimental import jitclass
from pandas.api.types import CategoricalDtype

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType


@jitclass
class RPacketLastInteractionTracker:
    r: nb.float64  # type: ignore[misc]
    nu: nb.float64  # type: ignore[misc]
    mu: nb.float64  # type: ignore[misc]
    energy: nb.float64  # type: ignore[misc]
    shell_id: nb.int64  # type: ignore[misc]
    interaction_type: nb.int64  # type: ignore[misc]

    # Interaction tracking (before/after)
    before_nu: nb.float64  # type: ignore[misc]
    before_mu: nb.float64  # type: ignore[misc]
    before_energy: nb.float64  # type: ignore[misc]
    after_nu: nb.float64  # type: ignore[misc]
    after_mu: nb.float64  # type: ignore[misc]
    after_energy: nb.float64  # type: ignore[misc]

    # Line interaction IDs (for line interactions only)
    interaction_line_absorb_id: nb.int64  # type: ignore[misc]
    interaction_line_emit_id: nb.int64  # type: ignore[misc]

    # Interaction counter
    interactions_count: nb.int64  # type: ignore[misc]

    """
    Numba JITCLASS for storing the last interaction the RPacket undergoes,
    with unified tracking for all interaction types.

    Parameters
    ----------
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
        interaction_before_* : float
            Properties before interaction (nu, mu, energy)
        interaction_after_* : float
            Properties after interaction (nu, mu, energy)
        interaction_line_absorb_id : int
            Line ID for absorbed line interactions (-1 for non-line interactions)
        interaction_line_emit_id : int
            Line ID for emitted line interactions (-1 for non-line interactions)
        interactions_count : int
            Count of interactions
    """

    def __init__(self) -> None:
        """
        Initialize properties with default values.
        Float values are initialized with NaN for better data quality tracking.
        """
        self.r = -1.0
        self.nu = float('nan')
        self.energy = float('nan')
        self.shell_id = -1
        self.interaction_type = -1

        # Initialize interaction tracking
        self.before_nu = float('nan')
        self.before_mu = float('nan')
        self.before_energy = float('nan')
        self.after_nu = float('nan')
        self.after_mu = float('nan')
        self.after_energy = float('nan')

        # Initialize line interaction IDs
        self.interaction_line_absorb_id = -1
        self.interaction_line_emit_id = -1

        # Initialize counter
        self.interactions_count = 0

    def track_line_interaction_before(self, r_packet) -> None:
        """
        Track packet state before line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before line interaction.
        """
        self.before_nu = r_packet.nu
        self.before_mu = r_packet.mu
        self.before_energy = r_packet.energy
        # Track the line ID that will be absorbed
        self.interaction_line_absorb_id = r_packet.next_line_id

    def track_line_interaction_after(self, r_packet) -> None:
        """
        Track packet state after line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after line interaction.
        """
        self.after_nu = r_packet.nu
        self.after_mu = r_packet.mu
        self.after_energy = r_packet.energy
        # Track the line ID that was emitted
        self.interaction_line_emit_id = r_packet.last_line_interaction_out_id
        self.interactions_count += 1
        # Update general tracking
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = r_packet.last_interaction_type

    def track_escattering_interaction_before(self, r_packet) -> None:
        """
        Track packet state before electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before electron scattering.
        """
        self.before_mu = r_packet.mu
        self.before_nu = r_packet.nu
        self.before_energy = r_packet.energy
        # Reset line IDs for non-line interactions
        self.interaction_line_absorb_id = -1
        self.interaction_line_emit_id = -1

    def track_escattering_interaction_after(self, r_packet) -> None:
        """
        Track packet state after electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after electron scattering.
        """
        self.after_mu = r_packet.mu
        self.after_nu = r_packet.nu
        self.after_energy = r_packet.energy
        self.interactions_count += 1
        # Update general tracking
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = r_packet.last_interaction_type

    def track_continuum_interaction_before(self, r_packet) -> None:
        """
        Track packet state before continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before continuum process.
        """
        self.before_nu = r_packet.nu
        self.before_energy = r_packet.energy
        self.before_mu = r_packet.mu
        # Reset line IDs for non-line interactions
        self.interaction_line_absorb_id = -1
        self.interaction_line_emit_id = -1

    def track_continuum_interaction_after(self, r_packet) -> None:
        """
        Track packet state after continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after continuum process.
        """
        self.after_nu = r_packet.nu
        self.after_energy = r_packet.energy
        self.after_mu = r_packet.mu
        self.interactions_count += 1
        # Update general tracking
        self.r = r_packet.r
        self.nu = r_packet.nu
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = r_packet.last_interaction_type

    def track_boundary_event(self, r_packet, from_shell_id=-1, to_shell_id=-1) -> None:
        """
        Track packet state during boundary event.
        This method provides API compatibility with RPacketTracker.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet during boundary event.
        from_shell_id : int, optional
            Shell ID the packet is leaving (default: -1).
        to_shell_id : int, optional
            Shell ID the packet is entering (default: -1).
        """
        pass
    
    def get_interaction_summary(self) -> nb.int64:  # type: ignore[misc]
        """
        Get summary of interaction count.

        Returns
        -------
        int
            Total interaction count.
        """
        return self.interactions_count


    def finalize_array(self) -> None:
        """
        Added to make RPacketLastInteractionTracker compatible with RPacketTracker
        """


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
    # Create individual trackers - the List will be created externally
    trackers = []
    for i in range(no_of_packets):
        trackers.append(RPacketLastInteractionTracker())
    return trackers


def rpacket_last_interaction_tracker_list_to_dataframe(tracker_list):
    """
    Convert a list of RPacketLastInteractionTracker instances to a DataFrame.

    This function extracts the last interaction data from each tracker and
    creates a pandas DataFrame with the relevant interaction information.

    Parameters
    ----------
    tracker_list : list
        List of RPacketLastInteractionTracker instances

    Returns
    -------
    pd.DataFrame
        DataFrame containing last interaction data with columns:
        - last_interaction_type: Type of the last interaction (categorical)
        - before_nu: Frequency before interaction
        - before_mu: Direction cosine before interaction
        - before_energy: Energy before interaction
        - after_nu: Frequency after interaction
        - after_mu: Direction cosine after interaction
        - after_energy: Energy after interaction
        - line_absorb_id: Line ID for absorbed line interactions (-1 for non-line)
        - line_emit_id: Line ID for emitted line interactions (-1 for non-line)
        - interactions_count: Total interaction count
    """
    # Create categorical dtype from enum for better performance and type safety
    interaction_type_dtype = CategoricalDtype(
        categories=[member.name for member in InteractionType],
        ordered=False
    )

    # Extract data from trackers
    interaction_type_raw = np.array([tracker.interaction_type for tracker in tracker_list])

    # Convert enum values to their string names and create categorical
    interaction_type_labels = [InteractionType(int_type).name for int_type in interaction_type_raw]
    last_interaction_type = pd.Categorical(
        interaction_type_labels,
        dtype=interaction_type_dtype
    )

    # Extract interaction data
    interaction_before_nu = np.array([tracker.interaction_before_nu for tracker in tracker_list])
    interaction_before_mu = np.array([tracker.interaction_before_mu for tracker in tracker_list])
    interaction_before_energy = np.array([tracker.interaction_before_energy for tracker in tracker_list])
    interaction_after_nu = np.array([tracker.interaction_after_nu for tracker in tracker_list])
    interaction_after_mu = np.array([tracker.interaction_after_mu for tracker in tracker_list])
    interaction_after_energy = np.array([tracker.interaction_after_energy for tracker in tracker_list])
    interaction_line_absorb_id = np.array([tracker.interaction_line_absorb_id for tracker in tracker_list])
    interaction_line_emit_id = np.array([tracker.interaction_line_emit_id for tracker in tracker_list])
    interactions_count = np.array([tracker.interactions_count for tracker in tracker_list])

    # Create DataFrame with packet index
    df = pd.DataFrame({
        "last_interaction_type": last_interaction_type,
        "before_nu": interaction_before_nu,
        "before_mu": interaction_before_mu,
        "before_energy": interaction_before_energy,
        "after_nu": interaction_after_nu,
        "after_mu": interaction_after_mu,
        "after_energy": interaction_after_energy,
        "line_absorb_id": interaction_line_absorb_id,
        "line_emit_id": interaction_line_emit_id,
        "interactions_count": interactions_count,
    }, index=pd.RangeIndex(len(tracker_list), name="packet_id"))

    return df
