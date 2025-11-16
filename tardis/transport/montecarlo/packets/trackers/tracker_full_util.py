import numpy as np
import pandas as pd
from numba import njit
from pandas.api.types import CategoricalDtype

from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    TrackerFull,
)


@njit
def trackers_full_list_to_arrays(trackers) -> tuple:
    """
    Convert a list of TrackerFull instances to concatenated arrays efficiently.

    Parameters
    ----------
    trackers : numba.typed.typedlist.List
        List of TrackerFull instances.

    Returns
    -------
    tuple
        Tuple containing:
        - packet_id (np.int64[:]) : Packet ID for each event
        - event_id (np.int64[:]) : Event ID within each packet
        - r (np.float64[:]) : Radius at each event
        - before_shell_id (np.int64[:]) : Shell ID before each event
        - after_shell_id (np.int64[:]) : Shell ID after each event
        - interaction_type (np.int64[:]) : Interaction type at each event
        - status (np.int64[:]) : Status at each event
        - before_nu (np.float64[:]) : Frequency before interaction
        - before_mu (np.float64[:]) : Mu before interaction
        - before_energy (np.float64[:]) : Energy before interaction
        - after_nu (np.float64[:]) : Frequency after interaction
        - after_mu (np.float64[:]) : Mu after interaction
        - after_energy (np.float64[:]) : Energy after interaction
        - line_absorb_id (np.int64[:]) : Line absorb ID
        - line_emit_id (np.int64[:]) : Line emit ID
    """
    # Calculate total array size
    total_events = 0
    for tracker in trackers:
        total_events += len(tracker.radius)

    # Pre-allocate all arrays
    packet_id = np.empty(total_events, dtype=np.int64)
    event_id = np.empty(total_events, dtype=np.int64)
    radius = np.empty(total_events, dtype=np.float64)
    before_shell_id = np.empty(total_events, dtype=np.int64)
    after_shell_id = np.empty(total_events, dtype=np.int64)
    interaction_type = np.empty(total_events, dtype=np.int64)
    status = np.empty(total_events, dtype=np.int64)
    before_nu = np.empty(total_events, dtype=np.float64)
    before_mu = np.empty(total_events, dtype=np.float64)
    before_energy = np.empty(total_events, dtype=np.float64)
    after_nu = np.empty(total_events, dtype=np.float64)
    after_mu = np.empty(total_events, dtype=np.float64)
    after_energy = np.empty(total_events, dtype=np.float64)
    line_absorb_id = np.empty(total_events, dtype=np.int64)
    line_emit_id = np.empty(total_events, dtype=np.int64)

    # Concatenate arrays from all trackers
    offset = 0
    for packet_idx, tracker in enumerate(trackers):
        n_events = len(tracker.radius)
        end_idx = offset + n_events

        # Simple array concatenation
        packet_id[offset:end_idx] = packet_idx
        event_id[offset:end_idx] = np.arange(n_events, dtype=np.int64)
        radius[offset:end_idx] = tracker.radius
        before_shell_id[offset:end_idx] = tracker.before_shell_id
        after_shell_id[offset:end_idx] = tracker.after_shell_id
        interaction_type[offset:end_idx] = tracker.interaction_type
        status[offset:end_idx] = tracker.status
        before_nu[offset:end_idx] = tracker.before_nu
        before_mu[offset:end_idx] = tracker.before_mu
        before_energy[offset:end_idx] = tracker.before_energy
        after_nu[offset:end_idx] = tracker.after_nu
        after_mu[offset:end_idx] = tracker.after_mu
        after_energy[offset:end_idx] = tracker.after_energy
        line_absorb_id[offset:end_idx] = tracker.line_absorb_id
        line_emit_id[offset:end_idx] = tracker.line_emit_id

        offset = end_idx

    return (
        packet_id,
        event_id,
        radius,
        before_shell_id,
        after_shell_id,
        interaction_type,
        status,
        before_nu,
        before_mu,
        before_energy,
        after_nu,
        after_mu,
        after_energy,
        line_absorb_id,
        line_emit_id,
    )


def trackers_full_to_df(trackers_full_list) -> pd.DataFrame:
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
    # Use the fast array extraction function
    (
        packet_id,
        event_id,
        r,
        before_shell_id,
        after_shell_id,
        interaction_type,
        status,
        before_nu,
        before_mu,
        before_energy,
        after_nu,
        after_mu,
        after_energy,
        line_absorb_id,
        line_emit_id,
    ) = trackers_full_list_to_arrays(trackers_full_list)

    # Create categorical data types
    interaction_type_dtype = CategoricalDtype(
        categories=[member.name for member in InteractionType], ordered=False
    )
    status_dtype = CategoricalDtype(
        categories=[member.name for member in PacketStatus], ordered=False
    )

    # Convert enum values to their string names and create categorical
    interaction_type_labels = [
        InteractionType(int(int_type)).name for int_type in interaction_type
    ]
    interaction_type_categorical = pd.Categorical(
        interaction_type_labels, dtype=interaction_type_dtype
    )

    # Convert status values to their string names
    status_labels = [
        PacketStatus(int(int_status)).name for int_status in status
    ]
    status_categorical = pd.Categorical(status_labels, dtype=status_dtype)

    # Create comprehensive event DataFrame with specific column order
    event_data = {
        "interaction_type": interaction_type_categorical,
        "status": status_categorical,
        "radius": r,
        "before_shell_id": before_shell_id,
        "after_shell_id": after_shell_id,
        "before_nu": before_nu,
        "before_mu": before_mu,
        "before_energy": before_energy,
        "after_nu": after_nu,
        "after_mu": after_mu,
        "after_energy": after_energy,
        "line_absorb_id": pd.array(line_absorb_id, dtype="int64"),
        "line_emit_id": pd.array(line_emit_id, dtype="int64"),
    }

    # Create multi-index with packet_id and event_id
    multi_index = pd.MultiIndex.from_arrays(
        [packet_id, event_id], names=["packet_id", "event_id"]
    )

    return pd.DataFrame(event_data, index=multi_index)


def tracker_full_df2tracker_last_interaction_df(
    full_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Convert a full tracking dataframe to show only the last interaction for each packet.

    Takes a full event tracking dataframe (from trackers_full_to_dataframe) and
    filters it to show only the last physics interaction each packet experienced.

    Parameters
    ----------
    full_df : pd.DataFrame
        Full event tracking dataframe with multi-index (packet_id, event_id) and
        columns ['radius', 'before_shell_id', 'after_shell_id', 'interaction_type'].

    Returns
    -------
    pd.DataFrame
        DataFrame containing last interaction data for each packet,
        indexed by packet_id with columns ['radius', 'shell_id', 'last_interaction_type',
        'line_absorb_id', 'line_emit_id', 'before_nu', 'before_mu', 'before_energy',
        'after_nu', 'after_mu', 'after_energy']. Note: shell_id comes from before_shell_id.
    """
    if full_df.empty:
        return pd.DataFrame()

    # Reset index to work with packet_id and event_id as columns
    df = full_df.reset_index()

    # Get all unique packet IDs to ensure we include packets with no physics interactions
    all_packet_ids = df["packet_id"].unique()

    # Filter out boundary events to only get actual interactions
    physics_interactions = df[
        df["interaction_type"].isin(
            ["LINE", "ESCATTERING", "CONTINUUM_PROCESS"]
        )
    ]

    if physics_interactions.empty:
        # No physics interactions - create DataFrame with NO_INTERACTION for all packets
        no_interaction_df = pd.DataFrame(
            {
                "last_interaction_type": pd.Categorical(
                    ["NO_INTERACTION"] * len(all_packet_ids)
                ),
                "status": pd.Categorical(["IN_PROCESS"] * len(all_packet_ids)),
                "radius": np.full(len(all_packet_ids), np.nan),
                "shell_id": np.full(len(all_packet_ids), -1),
                "before_nu": np.full(len(all_packet_ids), np.nan),
                "before_mu": np.full(len(all_packet_ids), np.nan),
                "before_energy": np.full(len(all_packet_ids), np.nan),
                "after_nu": np.full(len(all_packet_ids), np.nan),
                "after_mu": np.full(len(all_packet_ids), np.nan),
                "after_energy": np.full(len(all_packet_ids), np.nan),
                "line_absorb_id": pd.array(
                    np.full(len(all_packet_ids), -1), dtype="int64"
                ),
                "line_emit_id": pd.array(
                    np.full(len(all_packet_ids), -1), dtype="int64"
                ),
            },
            index=pd.Index(all_packet_ids, name="packet_id"),
        )
        return no_interaction_df

    # Get the last interaction for each packet
    last_interactions = physics_interactions.groupby("packet_id").last()

    # Rename columns to match last interaction format
    last_interaction_df = last_interactions.copy()
    last_interaction_df.rename(
        columns={
            "interaction_type": "last_interaction_type",
            "before_shell_id": "shell_id",
        },
        inplace=True,
    )

    # Drop after_shell_id column
    last_interaction_df = last_interaction_df.drop(columns=["after_shell_id"])

    # Handle packets that had no physics interactions - add them with NO_INTERACTION
    packets_with_interactions = set(last_interaction_df.index)
    packets_without_interactions = (
        set(all_packet_ids) - packets_with_interactions
    )

    if packets_without_interactions:
        # Create DataFrame for packets with no interactions
        no_interaction_data = {
            "last_interaction_type": pd.Categorical(
                ["NO_INTERACTION"] * len(packets_without_interactions)
            ),
            "status": pd.Categorical(
                ["IN_PROCESS"] * len(packets_without_interactions)
            ),
            "radius": np.full(len(packets_without_interactions), np.nan),
            "shell_id": np.full(len(packets_without_interactions), -1),
            "before_nu": np.full(len(packets_without_interactions), np.nan),
            "before_mu": np.full(len(packets_without_interactions), np.nan),
            "before_energy": np.full(len(packets_without_interactions), np.nan),
            "after_nu": np.full(len(packets_without_interactions), np.nan),
            "after_mu": np.full(len(packets_without_interactions), np.nan),
            "after_energy": np.full(len(packets_without_interactions), np.nan),
            "line_absorb_id": pd.array(
                np.full(len(packets_without_interactions), -1), dtype="int64"
            ),
            "line_emit_id": pd.array(
                np.full(len(packets_without_interactions), -1), dtype="int64"
            ),
        }
        no_interaction_df = pd.DataFrame(
            no_interaction_data,
            index=pd.Index(
                list(packets_without_interactions), name="packet_id"
            ),
        )

        # Combine the DataFrames
        last_interaction_df = pd.concat(
            [last_interaction_df, no_interaction_df], ignore_index=False
        )
        last_interaction_df.sort_index(inplace=True)

    return last_interaction_df


@njit
def generate_tracker_full_list(
    no_of_packets: int, length: int, extend_factor: int = 2
) -> list:
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
        rpacket_trackers.append(
            TrackerFull(length, extend_factor)
        )
    return rpacket_trackers
