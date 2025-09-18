import numpy as np
import pandas as pd
from numba import njit
from pandas.api.types import CategoricalDtype

from tardis.transport.montecarlo.packets.radiative_packet import (
    InteractionType,
    PacketStatus,
)
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import (
    TrackerLastInteraction,
)


@njit
def generate_tracker_last_interaction_list(no_of_packets):
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
        trackers.append(TrackerLastInteraction())
    return trackers


def trackers_last_interaction_to_df(tracker_list):
    """
    Convert a list of RPacketLastInteractionTracker instances to a DataFrame.

    This function extracts the last interaction data from each tracker and
    creates a pandas DataFrame with the same structure as tracker_full_df2tracker_last_interaction_df.

    Parameters
    ----------
    tracker_list : list
        List of RPacketLastInteractionTracker instances

    Returns
    -------
    pd.DataFrame
        DataFrame containing last interaction data with columns matching
        tracker_full_df2tracker_last_interaction_df output:
        - last_interaction_type: Type of the last interaction (categorical)
        - status: Packet status (categorical)
        - radius: Radius at interaction (NaN for last interaction trackers)
        - shell_id: Shell ID where interaction occurred (-1 for last interaction trackers)
        - before_nu: Frequency before interaction
        - before_mu: Direction cosine before interaction
        - before_energy: Energy before interaction
        - after_nu: Frequency after interaction
        - after_mu: Direction cosine after interaction
        - after_energy: Energy after interaction
                - line_absorb_id: Line ID for absorbed line interactions (-1 for non-line)
        - line_emit_id: Line ID for emitted line interactions (-1 for non-line)
    """
    # Create categorical dtypes from enums for better performance and type safety
    """
    from tardis.transport.montecarlo.packets.radiative_packet import PacketStatus
    
        """
    # Create categorical dtypes from enums for better performance and type safety
    interaction_type_dtype = CategoricalDtype(
        categories=[member.name for member in InteractionType], ordered=False
    )
    status_dtype = CategoricalDtype(
        categories=[member.name for member in PacketStatus], ordered=False
    )

    # Extract data from trackers
    interaction_type_raw = np.array(
        [tracker.interaction_type for tracker in tracker_list]
    )

    # Convert enum values to their string names and create categorical
    # Handle NO_INTERACTION (-1) specifically
    interaction_type_labels = []
    for int_type in interaction_type_raw:
        if int_type == -1:  # NO_INTERACTION
            interaction_type_labels.append("NO_INTERACTION")
        else:
            interaction_type_labels.append(InteractionType(int_type).name)
    
    last_interaction_type = pd.Categorical(
        interaction_type_labels, dtype=interaction_type_dtype
    )

    # Create status categorical (all IN_PROCESS for last interaction trackers)
    status_labels = ["IN_PROCESS"] * len(tracker_list)
    status_categorical = pd.Categorical(status_labels, dtype=status_dtype)

    # Extract interaction data
    before_nu = np.array([tracker.before_nu for tracker in tracker_list])
    before_mu = np.array([tracker.before_mu for tracker in tracker_list])
    before_energy = np.array([tracker.before_energy for tracker in tracker_list])
    after_nu = np.array([tracker.after_nu for tracker in tracker_list])
    after_mu = np.array([tracker.after_mu for tracker in tracker_list])
    after_energy = np.array([tracker.after_energy for tracker in tracker_list])
    line_absorb_id = np.array([tracker.interaction_line_absorb_id for tracker in tracker_list])
    line_emit_id = np.array([tracker.interaction_line_emit_id for tracker in tracker_list])

    # Create DataFrame with same column order as tracker_full_df2tracker_last_interaction_df
    df = pd.DataFrame(
        {
            "last_interaction_type": last_interaction_type,
            "status": status_categorical,
            "radius": np.full(len(tracker_list), np.nan),  # Not available in last interaction tracker
            "shell_id": np.full(len(tracker_list), -1),    # Not available in last interaction tracker
            "before_nu": before_nu,
            "before_mu": before_mu,
            "before_energy": before_energy,
            "after_nu": after_nu,
            "after_mu": after_mu,
            "after_energy": after_energy,
            "line_absorb_id": pd.array(line_absorb_id, dtype="int64"),
            "line_emit_id": pd.array(line_emit_id, dtype="int64"),
        },
        index=pd.RangeIndex(len(tracker_list), name="packet_id"),
    )

    return df