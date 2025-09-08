import numba as nb
import numpy as np
import pandas as pd
from numba import njit
from numba.experimental import jitclass
from pandas.api.types import CategoricalDtype

from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
)
from tardis.transport.montecarlo.packets.trackers.array_utils import (
    extend_array,
    extend_float_array,
    extend_int_array,
    extend_interaction_type_array,
)

NO_INTERACTION_INT = int(InteractionType.NO_INTERACTION)


@jitclass
class RPacketTracker:
    """
    Numba JITCLASS for storing interaction information for RPacket instances.

    Tracks core packet state (r, shell_id, interaction_type) for all events,
    plus detailed before/after interaction data with linked event IDs.

    Parameters
    ----------
    length : int
        Initial length of tracking arrays.

    Attributes
    ----------
    r : nb.float64[:]
        Radius of the RPacket at each event.
    shell_id : nb.int64[:]
        Shell ID where the RPacket is located at each event.
    interaction_type : nb.int64[:]
        Type of interaction undergone by the RPacket at each event.
    status : nb.int64[:]
        Status of the RPacket at each event.
    interaction_event_id : nb.int64[:]
        Event ID linked to each interaction (links before/after data).
    energy : nb.float64[:]
        Energy of the RPacket at each interaction.
    interaction_before_nu : nb.float64[:]
        Frequency before each interaction.
    interaction_before_mu : nb.float64[:]
        Cosine of angle before each interaction.
    interaction_after_nu : nb.float64[:]
        Frequency after each interaction.
    interaction_after_mu : nb.float64[:]
        Cosine of angle after each interaction.
    interaction_line_absorb_id : nb.int64[:]
        Line ID for absorbed line interactions.
    interaction_line_emit_id : nb.int64[:]
        Line ID for emitted line interactions.
    event_id : nb.int64
        Current event counter.
    current_interaction_index : nb.int64
        Current interaction counter.
    extend_factor : nb.int64
        Factor by which to extend arrays when capacity is reached.
    """

    # Core arrays that track everything (every event)
    r: nb.float64[:]  # type: ignore[misc]
    shell_id: nb.int64[:]  # type: ignore[misc]
    interaction_type: nb.int64[:]  # type: ignore[misc]
    status: nb.int64[:]  # type: ignore[misc]

    # Interaction-specific arrays with event ID linking
    interaction_event_id: nb.int64[:]  # type: ignore[misc]
    energy: nb.float64[:]  # type: ignore[misc]
    interaction_before_nu: nb.float64[:]  # type: ignore[misc]
    interaction_before_mu: nb.float64[:]  # type: ignore[misc]
    interaction_after_nu: nb.float64[:]  # type: ignore[misc]
    interaction_after_mu: nb.float64[:]  # type: ignore[misc]

    # Line interaction specific tracking
    interaction_line_absorb_id: nb.int64[:]  # type: ignore[misc]
    interaction_line_emit_id: nb.int64[:]  # type: ignore[misc]

    event_id: nb.int64  # type: ignore[misc]
    current_interaction_index: nb.int64  # type: ignore[misc]
    extend_factor: nb.int64  # type: ignore[misc]

    def __init__(self, length: int) -> None:
        """
        Initialize the RPacketTracker with arrays for tracking packet events.

        Parameters
        ----------
        length : int
            Initial length of the tracking arrays.
        """
        # Core arrays that track every event
        self.r = np.empty(length, dtype=np.float64)
        self.shell_id = np.empty(length, dtype=np.int64)
        self.interaction_type = np.full(
            length, NO_INTERACTION_INT, dtype=np.int64
        )
        self.status = np.empty(length, dtype=np.int64)

        # Interaction arrays
        self.interaction_event_id = np.empty(length, dtype=np.int64)
        self.energy = np.empty(length, dtype=np.float64)
        self.interaction_before_nu = np.full(length, np.nan, dtype=np.float64)
        self.interaction_before_mu = np.full(length, np.nan, dtype=np.float64)
        self.interaction_after_nu = np.full(length, np.nan, dtype=np.float64)
        self.interaction_after_mu = np.full(length, np.nan, dtype=np.float64)

        # Line interaction specific tracking
        self.interaction_line_absorb_id = np.full(length, -1, dtype=np.int64)
        self.interaction_line_emit_id = np.full(length, -1, dtype=np.int64)

        self.event_id = 0
        self.current_interaction_index = 0
        self.extend_factor = 2

    @property
    def array_length(self) -> int:
        """
        Current capacity of core tracking arrays.

        Returns
        -------
        int
            Current array capacity.
        """
        return self.r.size

    @property
    def interaction_array_length(self) -> int:
        """
        Current capacity of interaction tracking arrays.

        Returns
        -------
        int
            Current interaction array capacity.
        """
        return self.interaction_event_id.size

    def _extend_core_arrays(self) -> None:
        """Extend core tracking arrays (r, shell_id, interaction_type, status)."""
        new_length = self.array_length * self.extend_factor

        self.r = extend_array(self.r, new_length)
        self.shell_id = extend_array(self.shell_id, new_length)
        self.interaction_type = extend_interaction_type_array(
            self.interaction_type, new_length
        )
        self.status = extend_array(self.status, new_length)

    def _extend_interaction_arrays(self) -> None:
        """Extend interaction tracking arrays."""
        new_length = self.interaction_array_length * self.extend_factor

        self.interaction_event_id = extend_array(
            self.interaction_event_id, new_length
        )
        self.energy = extend_array(self.energy, new_length)
        self.interaction_before_nu = extend_float_array(
            self.interaction_before_nu, new_length
        )
        self.interaction_before_mu = extend_float_array(
            self.interaction_before_mu, new_length
        )
        self.interaction_after_nu = extend_float_array(
            self.interaction_after_nu, new_length
        )
        self.interaction_after_mu = extend_float_array(
            self.interaction_after_mu, new_length
        )

        # Line interaction specific tracking
        self.interaction_line_absorb_id = extend_int_array(
            self.interaction_line_absorb_id, new_length
        )
        self.interaction_line_emit_id = extend_int_array(
            self.interaction_line_emit_id, new_length
        )

    def track_line_interaction_before(self, r_packet) -> None:
        """
        Track packet state before line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before line interaction.
        """
        # Extend core arrays if needed
        if self.event_id >= self.array_length:
            self._extend_core_arrays()

        # Extend interaction arrays if needed
        if self.current_interaction_index >= self.interaction_array_length:
            self._extend_interaction_arrays()

        # Track general packet state
        self.r[self.event_id] = r_packet.r
        self.shell_id[self.event_id] = r_packet.current_shell_id
        self.interaction_type[self.event_id] = r_packet.last_interaction_type
        self.status[self.event_id] = r_packet.status

        # Track interaction data
        self.interaction_event_id[self.current_interaction_index] = self.event_id
        self.energy[self.current_interaction_index] = r_packet.energy
        self.interaction_before_nu[self.current_interaction_index] = r_packet.nu
        self.interaction_before_mu[self.current_interaction_index] = r_packet.mu
        self.interaction_line_absorb_id[self.current_interaction_index] = (
            r_packet.next_line_id
        )

    def track_line_interaction_after(self, r_packet) -> None:
        """
        Track packet state after line interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after line interaction.
        """
        # Track line interaction after state for the SAME interaction_id
        self.interaction_after_nu[self.current_interaction_index] = r_packet.nu
        self.interaction_after_mu[self.current_interaction_index] = r_packet.mu
        self.interaction_line_emit_id[self.current_interaction_index] = (
            r_packet.last_line_interaction_out_id
        )

        # Increment counters
        self.event_id += 1
        self.current_interaction_index += 1

    def track_escattering_interaction_before(self, r_packet) -> None:
        """
        Track packet state before electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before electron scattering.
        """
        # Extend core arrays if needed
        if self.event_id >= self.array_length:
            self._extend_core_arrays()

        # Extend interaction arrays if needed
        if self.current_interaction_index >= self.interaction_array_length:
            self._extend_interaction_arrays()

        # Track general packet state
        self.r[self.event_id] = r_packet.r
        self.shell_id[self.event_id] = r_packet.current_shell_id
        self.interaction_type[self.event_id] = r_packet.last_interaction_type
        self.status[self.event_id] = r_packet.status

        # Track interaction data
        self.interaction_event_id[self.current_interaction_index] = self.event_id
        self.energy[self.current_interaction_index] = r_packet.energy
        self.interaction_before_nu[self.current_interaction_index] = r_packet.nu
        self.interaction_before_mu[self.current_interaction_index] = r_packet.mu

    def track_escattering_interaction_after(self, r_packet) -> None:
        """
        Track packet state after electron scattering interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after electron scattering.
        """
        # Track electron scattering after state for the SAME interaction_id
        # For Thomson scattering, frequency doesn't change
        self.interaction_after_nu[self.current_interaction_index] = r_packet.nu
        self.interaction_after_mu[self.current_interaction_index] = r_packet.mu

        # Increment counters
        self.event_id += 1
        self.current_interaction_index += 1

    def track_continuum_interaction_before(self, r_packet) -> None:
        """
        Track packet state before continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet before continuum process.
        """
        # Extend core arrays if needed
        if self.event_id >= self.array_length:
            self._extend_core_arrays()

        # Extend interaction arrays if needed
        if self.current_interaction_index >= self.interaction_array_length:
            self._extend_interaction_arrays()

        # Track general packet state
        self.r[self.event_id] = r_packet.r
        self.shell_id[self.event_id] = r_packet.current_shell_id
        self.interaction_type[self.event_id] = r_packet.last_interaction_type
        self.status[self.event_id] = r_packet.status

        # Track interaction data
        self.interaction_event_id[self.current_interaction_index] = self.event_id
        self.energy[self.current_interaction_index] = r_packet.energy
        self.interaction_before_nu[self.current_interaction_index] = r_packet.nu
        self.interaction_before_mu[self.current_interaction_index] = r_packet.mu

    def track_continuum_interaction_after(self, r_packet) -> None:
        """
        Track packet state after continuum process interaction.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet after continuum process.
        """
        # Track continuum after state for the SAME interaction_id
        self.interaction_after_nu[self.current_interaction_index] = r_packet.nu
        self.interaction_after_mu[self.current_interaction_index] = r_packet.mu

        # Increment counters
        self.event_id += 1
        self.current_interaction_index += 1

    def track_boundary_event(
        self, r_packet, from_shell_id: int, to_shell_id: int
    ) -> None:
        """
        Track boundary event when packet crosses shell boundary.

        Parameters
        ----------
        r_packet : RPacket
            The R-packet crossing the boundary.
        from_shell_id : int
            Shell ID the packet is leaving.
        to_shell_id : int
            Shell ID the packet is entering.
        """
        # Extend core arrays if needed
        if self.event_id >= self.array_length:
            self._extend_core_arrays()

        # Track boundary event state
        self.r[self.event_id] = r_packet.r
        self.shell_id[self.event_id] = to_shell_id  # Record destination shell
        self.interaction_type[self.event_id] = InteractionType.BOUNDARY
        self.status[self.event_id] = r_packet.status

        self.event_id += 1

    def finalize_array(self) -> None:
        """
        Change the size of the arrays from length ( or multiple of length ) to
        the actual number of events/interactions.
        """
        # Finalize core arrays to actual event count
        self.r = self.r[: self.event_id]
        self.shell_id = self.shell_id[: self.event_id]
        self.interaction_type = self.interaction_type[: self.event_id]
        self.status = self.status[: self.event_id]

        # Finalize interaction arrays to actual interaction count
        self.interaction_event_id = self.interaction_event_id[
            : self.current_interaction_index
        ]
        self.energy = self.energy[: self.current_interaction_index]
        self.interaction_before_nu = self.interaction_before_nu[
            : self.current_interaction_index
        ]
        self.interaction_before_mu = self.interaction_before_mu[
            : self.current_interaction_index
        ]
        self.interaction_after_nu = self.interaction_after_nu[
            : self.current_interaction_index
        ]
        self.interaction_after_mu = self.interaction_after_mu[
            : self.current_interaction_index
        ]

        # Finalize line interaction specific tracking arrays
        self.interaction_line_absorb_id = self.interaction_line_absorb_id[
            : self.current_interaction_index
        ]
        self.interaction_line_emit_id = self.interaction_line_emit_id[
            : self.current_interaction_index
        ]


def rpacket_trackers_to_dataframe(rpacket_trackers) -> pd.DataFrame:
    """
    Convert a list of RPacketTracker instances to a comprehensive DataFrame.

    Parameters
    ----------
    rpacket_trackers : numba.typed.typedlist.List
        List of RPacketTracker instances.

    Returns
    -------
    pd.DataFrame
        DataFrame containing all event data including line interaction details,
        indexed by packet_id and event_id. Line interaction columns have default
        values (-1 for IDs, NaN for physical quantities) for non-interaction events.
    """
    interaction_type_dtype = CategoricalDtype(
        categories=[member.name for member in InteractionType], ordered=False
    )

    # Collect all core event data
    all_r = []
    all_shell_id = []
    all_interaction_type = []
    all_status = []
    all_packet_index = []
    all_event_index = []
    
    # Collect line interaction data
    all_line_absorb_id = []
    all_line_emit_id = []
    all_before_nu = []
    all_before_mu = []
    all_before_energy = []
    all_after_nu = []
    all_after_mu = []
    all_after_energy = []

    for i, rpacket_tracker in enumerate(rpacket_trackers):
        n_events = rpacket_tracker.event_id
        n_interactions = rpacket_tracker.current_interaction_index

        # Core event data (every event)
        all_r.append(rpacket_tracker.r[:n_events])
        all_shell_id.append(rpacket_tracker.shell_id[:n_events])
        all_interaction_type.append(rpacket_tracker.interaction_type[:n_events])
        all_status.append(rpacket_tracker.status[:n_events])
        all_packet_index.append(np.full(n_events, i))
        all_event_index.append(np.arange(n_events))
        
        # Initialize line interaction arrays with default values
        event_line_absorb_id = np.full(n_events, -1, dtype=np.int64)
        event_line_emit_id = np.full(n_events, -1, dtype=np.int64)
        event_before_nu = np.full(n_events, np.nan, dtype=np.float64)
        event_before_mu = np.full(n_events, np.nan, dtype=np.float64)
        event_before_energy = np.full(n_events, np.nan, dtype=np.float64)
        event_after_nu = np.full(n_events, np.nan, dtype=np.float64)
        event_after_mu = np.full(n_events, np.nan, dtype=np.float64)
        event_after_energy = np.full(n_events, np.nan, dtype=np.float64)
        
        # Fill in actual interaction data where available
        if n_interactions > 0:
            # Map interaction data to corresponding events using interaction_event_id
            for j in range(n_interactions):
                event_id = rpacket_tracker.interaction_event_id[j]
                if event_id < n_events:  # Safety check
                    event_line_absorb_id[event_id] = rpacket_tracker.interaction_line_absorb_id[j]
                    event_line_emit_id[event_id] = rpacket_tracker.interaction_line_emit_id[j]
                    event_before_nu[event_id] = rpacket_tracker.interaction_before_nu[j]
                    event_before_mu[event_id] = rpacket_tracker.interaction_before_mu[j]
                    event_before_energy[event_id] = rpacket_tracker.energy[j]
                    event_after_nu[event_id] = rpacket_tracker.interaction_after_nu[j]
                    event_after_mu[event_id] = rpacket_tracker.interaction_after_mu[j]
                    event_after_energy[event_id] = rpacket_tracker.energy[j]  # Energy same before/after for now
        
        # Append to collections
        all_line_absorb_id.append(event_line_absorb_id)
        all_line_emit_id.append(event_line_emit_id)
        all_before_nu.append(event_before_nu)
        all_before_mu.append(event_before_mu)
        all_before_energy.append(event_before_energy)
        all_after_nu.append(event_after_nu)
        all_after_mu.append(event_after_mu)
        all_after_energy.append(event_after_energy)

    # Concatenate core arrays
    combined_interaction_type = np.concatenate(all_interaction_type)
    combined_status = np.concatenate(all_status)

    # Convert enum values to their string names and create categorical
    interaction_type_labels = [
        InteractionType(int_type).name for int_type in combined_interaction_type
    ]
    interaction_type_categorical = pd.Categorical(
        interaction_type_labels, dtype=interaction_type_dtype
    )

    # Convert status values to their string names
    status_dtype = CategoricalDtype(
        categories=[member.name for member in PacketStatus], ordered=False
    )
    status_labels = [
        PacketStatus(int_status).name for int_status in combined_status
    ]
    status_categorical = pd.Categorical(status_labels, dtype=status_dtype)

    # Create comprehensive event DataFrame
    event_data = {
        "r": np.concatenate(all_r),
        "shell_id": np.concatenate(all_shell_id),
        "interaction_type": interaction_type_categorical,
        "status": status_categorical,
        "line_absorb_id": np.concatenate(all_line_absorb_id),
        "line_emit_id": np.concatenate(all_line_emit_id),
        "before_nu": np.concatenate(all_before_nu),
        "before_mu": np.concatenate(all_before_mu),
        "before_energy": np.concatenate(all_before_energy),
        "after_nu": np.concatenate(all_after_nu),
        "after_mu": np.concatenate(all_after_mu),
        "after_energy": np.concatenate(all_after_energy),
    }

    # Create multi-index with packet_id and event_id
    index_arrays = [
        np.concatenate(all_packet_index),
        np.concatenate(all_event_index),
    ]
    multi_index = pd.MultiIndex.from_arrays(
        index_arrays, names=["packet_id", "event_id"]
    )

    return pd.DataFrame(event_data, index=multi_index)


def rpacket_trackers_to_last_interaction_dataframe(
    rpacket_trackers,
) -> pd.DataFrame:
    """
    Convert RPacketTracker instances to a last interaction DataFrame.

    Extracts only the last interaction each packet experienced, filtering out
    boundary events to show only physics interactions (line, escattering, continuum).

    Parameters
    ----------
    rpacket_trackers : numba.typed.typedlist.List
        List of RPacketTracker instances.

    Returns
    -------
    pd.DataFrame
        DataFrame containing last interaction data for each packet,
        indexed by packet_id.
    """
    # First get the full event dataframe
    full_df = rpacket_trackers_to_dataframe(rpacket_trackers)

    if full_df.empty:
        return pd.DataFrame()

    # Reset index to work with packet_id and event_id as columns
    df = full_df.reset_index()

    # Filter out boundary events to only get actual interactions
    physics_interactions = df[
        df["interaction_type"].isin(
            ["LINE", "ESCATTERING", "CONTINUUM_PROCESS"]
        )
    ]

    if physics_interactions.empty:
        return pd.DataFrame()

    # Get the last interaction for each packet
    last_interactions = physics_interactions.groupby("packet_id").last()

    # Keep only the relevant columns and rename to match last interaction format
    last_interaction_df = last_interactions[
        ["r", "shell_id", "interaction_type"]
    ].copy()
    last_interaction_df.rename(
        columns={
            "interaction_type": "last_interaction_type",
            "shell_id": "last_shell_id",
        },
        inplace=True,
    )

    return last_interaction_df


def full_tracking_to_last_interaction_dataframe(
    full_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Convert a full tracking dataframe to show only the last interaction for each packet.

    Takes a full event tracking dataframe (from rpacket_trackers_to_dataframe) and
    filters it to show only the last physics interaction each packet experienced.

    Parameters
    ----------
    full_df : pd.DataFrame
        Full event tracking dataframe with multi-index (packet_id, event_id) and
        columns ['r', 'shell_id', 'interaction_type'].

    Returns
    -------
    pd.DataFrame
        DataFrame containing last interaction data for each packet,
        indexed by packet_id with columns ['r', 'last_shell_id', 'last_interaction_type'].
    """
    if full_df.empty:
        return pd.DataFrame()

    # Get all unique packet IDs from the full dataframe
    all_packet_ids = full_df.index.get_level_values(0).unique()

    # Remove all boundary events - filter for actual physics interactions only
    physics_interactions = full_df[
        full_df["interaction_type"].isin(
            ["LINE", "ESCATTERING", "CONTINUUM_PROCESS"]
        )
    ]

    if physics_interactions.empty:
        # No physics interactions - create DataFrame with NO_INTERACTION for all packets
        last_interaction_df = pd.DataFrame(
            {
                "r": np.nan,
                "last_shell_id": -1,
                "last_interaction_type": "NO_INTERACTION",
            },
            index=pd.Index(all_packet_ids, name="packet_id"),
        )
        return last_interaction_df

    # Group by packet_id and take the last interaction for each packet
    last_interactions = physics_interactions.groupby(level=0).last()

    # Rename columns to match last interaction format
    last_interaction_df = last_interactions.copy()
    last_interaction_df.rename(
        columns={
            "interaction_type": "last_interaction_type",
            "shell_id": "last_shell_id",
        },
        inplace=True,
    )

    # Handle packets that had no physics interactions
    missing_packets = all_packet_ids.difference(last_interaction_df.index)
    if len(missing_packets) > 0:
        no_interaction_data = pd.DataFrame(
            {
                "r": np.nan,
                "last_shell_id": -1,
                "last_interaction_type": "NO_INTERACTION",
            },
            index=pd.Index(missing_packets, name="packet_id"),
        )
        last_interaction_df = pd.concat(
            [last_interaction_df, no_interaction_data]
        )
        last_interaction_df.sort_index(inplace=True)

    return last_interaction_df


@njit
def generate_rpacket_tracker_list(no_of_packets: int, length: int):
    """
    Generate list of RPacketTracker instances.

    Parameters
    ----------
    no_of_packets : int
        The count of RPackets that are sent in the ejecta.
    length : int
        Initial length of the tracking array.

    Returns
    -------
    list
        A list containing RPacketTracker for each RPacket.
    """
    rpacket_trackers = []
    for i in range(no_of_packets):
        rpacket_trackers.append(RPacketTracker(length))
    return rpacket_trackers
