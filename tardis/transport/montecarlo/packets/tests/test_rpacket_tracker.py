import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType
from tardis.transport.montecarlo.packets.trackers.array_utils import (
    extend_array,
)

NO_INTERACTION_INT = int(InteractionType.NO_INTERACTION)


@pytest.fixture(scope="module")
def tracker_full_df(simulation_rpacket_tracking):
    "RPacketTracker object from the simulation"
    return simulation_rpacket_tracking.transport.transport_state.tracker_full_df


def test_extend_array():
    original_array = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    new_length = 10

    new_array = extend_array(original_array, new_length)

    assert new_array.size == new_length
    assert new_array.dtype == original_array.dtype
    npt.assert_allclose(original_array, new_array[:original_array.size])


def test_tracker_last_interaction_df_interaction_types(simulation_rpacket_tracking, regression_data):
    df = simulation_rpacket_tracking.transport.transport_state.tracker_last_interaction_df
    interaction_type_labels = df['last_interaction_type'].values
    interaction_types = np.array([InteractionType[label].value if label != 'NO_INTERACTION' else -1 for label in interaction_type_labels], dtype=np.int64)

    expected_types = regression_data.sync_ndarray(interaction_types)
    npt.assert_allclose(interaction_types, expected_types)


def test_tracker_last_interaction_df_shell_ids(simulation_rpacket_tracking, regression_data):
    df = simulation_rpacket_tracking.transport.transport_state.tracker_last_interaction_df
    interaction_type_labels = df['last_interaction_type'].values
    interaction_types = np.array([InteractionType[label].value if label != 'NO_INTERACTION' else -1 for label in interaction_type_labels], dtype=np.int64)
    mask = interaction_types == InteractionType.LINE
    shell_ids = df['shell_id'].values[mask]

    expected_shells = regression_data.sync_ndarray(shell_ids)
    npt.assert_allclose(shell_ids, expected_shells)


def test_tracker_last_interaction_df_output_nus(simulation_rpacket_tracking, regression_data):
    packet_collection = simulation_rpacket_tracking.transport.transport_state.packet_collection
    output_nus = packet_collection.output_nus

    expected_nus = regression_data.sync_ndarray(output_nus)
    npt.assert_allclose(output_nus, expected_nus)


def test_boundary_interactions(tracker_full_df, regression_data):
    """
    Validate boundary events per packet using the tracker_full_df.

    Extract per-packet event indices where interaction_type equals BOUNDARY
    from the tracker_full_df.
    """
    df = tracker_full_df

    packet_ids = df.index.get_level_values('packet_id').unique()
    no_of_packets = len(packet_ids)

    per_packet_event_ids = []
    max_len = 0
    for packet_id in packet_ids:
        packet_df = df.loc[packet_id]
        boundary_mask = packet_df['interaction_type'] == 'BOUNDARY'
        boundary_idx = packet_df[boundary_mask].index.values
        per_packet_event_ids.append(boundary_idx.astype(np.int64))
        if boundary_idx.size > max_len:
            max_len = boundary_idx.size

    if max_len == 0:
        max_len = 1
    obtained_boundary_interaction = np.full(
        (no_of_packets, max_len),
        -1,
        dtype=np.int64,
    )
    for i, event_ids in enumerate(per_packet_event_ids):
        obtained_boundary_interaction[i, : event_ids.size] = event_ids

    expected_boundary_interaction = regression_data.sync_ndarray(
        obtained_boundary_interaction
    )

    npt.assert_array_equal(obtained_boundary_interaction, expected_boundary_interaction)


def test_tracker_full_df_contents(tracker_full_df, regression_data):
    tracker_df = tracker_full_df.copy()
    tracker_df['interaction_type'] = tracker_df['interaction_type'].astype(str)
    tracker_df['status'] = tracker_df['status'].astype(str)

    expected = regression_data.sync_dataframe(tracker_df, key="tracker_full_df")
    pd.testing.assert_frame_equal(tracker_df, expected)


def test_tracker_last_interaction_df_contents(simulation_rpacket_tracking, regression_data):
    tracker_df = simulation_rpacket_tracking.transport.transport_state.tracker_last_interaction_df.copy()
    tracker_df['last_interaction_type'] = tracker_df['last_interaction_type'].astype(str)

    expected = regression_data.sync_dataframe(tracker_df, key="tracker_last_interaction_df")
    pd.testing.assert_frame_equal(tracker_df, expected)
