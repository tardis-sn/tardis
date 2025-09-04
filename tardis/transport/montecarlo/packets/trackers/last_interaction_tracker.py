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
    energy: nb.float64  # type: ignore[misc]
    shell_id: nb.int64  # type: ignore[misc]
    interaction_type: nb.int64  # type: ignore[misc]

    # LINE interaction tracking (before/after)
    line_before_nu: nb.float64  # type: ignore[misc]
    line_before_mu: nb.float64  # type: ignore[misc]
    line_before_id: nb.int64  # type: ignore[misc]
    line_after_nu: nb.float64  # type: ignore[misc]
    line_after_mu: nb.float64  # type: ignore[misc]
    line_after_id: nb.int64  # type: ignore[misc]

    # ESCATTERING interaction tracking (before/after)
    escat_before_mu: nb.float64  # type: ignore[misc]
    escat_after_mu: nb.float64  # type: ignore[misc]

    # CONTINUUM_PROCESS interaction tracking (before/after)
    continuum_before_nu: nb.float64  # type: ignore[misc]
    continuum_before_energy: nb.float64  # type: ignore[misc]
    continuum_before_mu: nb.float64  # type: ignore[misc]
    continuum_after_nu: nb.float64  # type: ignore[misc]
    continuum_after_energy: nb.float64  # type: ignore[misc]
    continuum_after_mu: nb.float64  # type: ignore[misc]

    # BOUNDARY interaction tracking (single state)
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
        line_before_* : float/int
            Properties before line interaction (nu, mu, id)
        line_after_* : float/int
            Properties after line interaction (nu, mu, id)
        escat_before_* : float
            Properties before electron scattering interaction (mu)
        escat_after_* : float
            Properties after electron scattering interaction (mu)
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
        Initialize properties with default values.
        Float values are initialized with NaN for better data quality tracking.
        """
        self.r = -1.0
        self.nu = float('nan')
        self.energy = float('nan')
        self.shell_id = -1
        self.interaction_type = -1

        # Initialize LINE interaction tracking
        self.line_before_nu = float('nan')
        self.line_before_mu = float('nan')
        self.line_before_id = -1
        self.line_after_nu = float('nan')
        self.line_after_mu = float('nan')
        self.line_after_id = -1

        # Initialize ESCATTERING interaction tracking
        self.escat_before_mu = float('nan')
        self.escat_after_mu = float('nan')

        # Initialize CONTINUUM_PROCESS interaction tracking
        self.continuum_before_nu = float('nan')
        self.continuum_before_energy = float('nan')
        self.continuum_before_mu = float('nan')
        self.continuum_after_nu = float('nan')
        self.continuum_after_energy = float('nan')
        self.continuum_after_mu = float('nan')

        # Initialize BOUNDARY interaction tracking
        self.boundary_r = -1.0
        self.boundary_nu = float('nan')
        self.boundary_energy = float('nan')
        self.boundary_mu = float('nan')
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
        self.line_before_mu = r_packet.mu
        # The line ID should already be set in last_line_interaction_in_id by the transport phase
        self.line_before_id = r_packet.last_line_interaction_in_id

    def track_line_interaction_after(self, r_packet) -> None:
        """
        Track packet state after line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after line interaction.
        """
        self.line_after_nu = r_packet.nu
        self.line_after_mu = r_packet.mu
        # Capture the line ID after the interaction (from last_line_interaction_out_id)
        self.line_after_id = r_packet.last_line_interaction_out_id
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
        self.escat_before_mu = r_packet.mu

    def track_escattering_interaction_after(self, r_packet) -> None:
        """
        Track packet state after electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after electron scattering.
        """
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
        self.boundary_mu = r_packet.mu
        self.boundary_from_shell_id = from_shell_id
        self.boundary_to_shell_id = to_shell_id
        self.boundary_interactions_count += 1
        # Update general tracking
        self.r = r_packet.r
        self.shell_id = r_packet.current_shell_id
        self.interaction_type = r_packet.last_interaction_type

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

    def initialize_tracker(self, r_packet):
        """
        Initialize tracker with packet properties.
        
        Parameters
        ----------
        r_packet : RPacket
            The R-packet to initialize tracking with
        """
        self.nu = r_packet.nu
        self.mu = r_packet.mu
        self.r = r_packet.r
        self.energy = r_packet.energy
        self.shell_id = r_packet.current_shell_id

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
        - last_line_interaction_in_id: Line ID for the input line interaction
        - last_line_interaction_in_nu: Frequency for the input line interaction
        - last_line_interaction_out_id: Line ID for the output line interaction
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

    # For line interactions, use line_before_nu as the input frequency
    # For other interactions, use the packet frequency
    last_line_interaction_in_nu = np.array([
        tracker.line_before_nu if tracker.line_interactions_count > 0 else tracker.nu
        for tracker in tracker_list
    ])

    # Extract line interaction IDs - only meaningful for LINE interactions
    last_line_interaction_in_id = np.array([
        tracker.line_before_id if tracker.interaction_type == InteractionType.LINE else -1
        for tracker in tracker_list
    ])
    last_line_interaction_out_id = np.array([
        tracker.line_after_id if tracker.interaction_type == InteractionType.LINE else -1
        for tracker in tracker_list
    ])

    # Create DataFrame with packet index
    df = pd.DataFrame({
        "last_interaction_type": last_interaction_type,
        "last_line_interaction_in_id": last_line_interaction_in_id,
        "last_line_interaction_in_nu": last_line_interaction_in_nu,
        "last_line_interaction_out_id": last_line_interaction_out_id,
    }, index=pd.RangeIndex(len(tracker_list), name="packet_id"))

    return df