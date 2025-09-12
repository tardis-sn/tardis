import numpy as np
import numpy.testing as npt
import pandas as pd
import pytest

from pathlib import Path
from tardis.transport.montecarlo.packets.radiative_packet import InteractionType
from tardis.transport.montecarlo.packets.trackers.array_utils import (
    extend_array,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full import (
    trackers_full_to_dataframe,
)

NO_INTERACTION_INT = int(InteractionType.NO_INTERACTION)


@pytest.fixture
def interaction_type_last_interaction_class(
    simulation_rpacket_tracking,
):
    """Last interaction types of rpacket from LastInteractionTracker class"""
    interaction_type = simulation_rpacket_tracking.transport.transport_state.last_interaction_type
    return interaction_type


@pytest.fixture
def shell_id_last_interaction_class(
    simulation_rpacket_tracking,
):
    """
    shell_id when last interaction is line from LastInteractionTracker class
    """
    interaction_type = simulation_rpacket_tracking.transport.transport_state.last_interaction_type
    mask = interaction_type == InteractionType.LINE
    shell_id = simulation_rpacket_tracking.transport.transport_state.last_line_interaction_shell_id
    last_line_interaction_shell_id = shell_id[mask]

    return last_line_interaction_shell_id


@pytest.fixture
def nu_from_packet_collection(
    simulation_rpacket_tracking,
):
    """Last interaction output nus of rpacket from packet_collection"""
    packet_collection = (
        simulation_rpacket_tracking.transport.transport_state.packet_collection
    )
    return packet_collection.output_nus


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


@pytest.mark.parametrize(
    "expected,obtained",
    [
        (
            "interaction_type_last_interaction_class",
            "last_interaction_type_rpacket_tracker",
        ),
        ("shell_id_last_interaction_class", "shell_id_rpacket_tracker"),
        ("nu_from_packet_collection", "nu_rpacket_tracker"),
    ],
)
def test_tracker_full_list_properties(expected, obtained, request):
    expected = request.getfixturevalue(expected)
    obtained = request.getfixturevalue(obtained)
    npt.assert_allclose(expected, obtained)


def test_boundary_interactions(tracker_full_df, regression_data):
    no_of_packets = len(tracker_full_df)
    expected_boundary_interaction = regression_data.sync_ndarray(
        obtained_boundary_interaction
    )

    max_boundary_interaction_size = max(
        [tracker.boundary_interaction.size for tracker in tracker_full_df]
    )
    obtained_boundary_interaction = np.full(
        (no_of_packets, max_boundary_interaction_size),
        [-1],
        dtype=tracker_full_df[0].boundary_interaction.dtype,
    )

    for i, tracker in enumerate(tracker_full_df):
        obtained_boundary_interaction[
            i, : tracker.boundary_interaction.size
        ] = tracker.boundary_interaction

    
    npt.assert_array_equal(
        obtained_boundary_interaction, expected_boundary_interaction
    )


def test_tracker_full_lists_to_dataframe(simulation_rpacket_tracking):
    transport_state = simulation_rpacket_tracking.transport.transport_state
    rtracker_df = trackers_full_to_dataframe(transport_state.rpacket_tracker)

    # check df shape and column names
    assert rtracker_df.shape == (
        sum([len(tracker.r) for tracker in transport_state.rpacket_tracker]),
        8,
    )
    npt.assert_array_equal(
        transport_state.rpacket_tracker_df.columns.values,
        np.array(
            [
                "status",
                "seed",
                "r",
                "nu",
                "mu",
                "energy",
                "shell_id",
                "interaction_type",
            ]
        ),
    )

    # check all data with rpacket_tracker
    expected_rtrackers = []
    for rpacket in transport_state.rpacket_tracker:
        for rpacket_step_no in range(len(rpacket.r)):
            expected_rtrackers.append(
                [
                    rpacket.status[rpacket_step_no],
                    rpacket.seed,
                    rpacket.r[rpacket_step_no],
                    rpacket.nu[rpacket_step_no],
                    rpacket.mu[rpacket_step_no],
                    rpacket.energy[rpacket_step_no],
                    rpacket.shell_id[rpacket_step_no],
                    rpacket.interaction_type[rpacket_step_no],
                ]
            )
    npt.assert_array_equal(rtracker_df.to_numpy(), np.array(expected_rtrackers))
