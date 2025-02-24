import numpy as np
import numpy.testing as npt
import pytest

from tardis.transport.montecarlo.packet_trackers import (
    RPacketTracker,
    rpacket_trackers_to_dataframe,
)
from tardis.transport.montecarlo.r_packet import InteractionType


@pytest.fixture
def interaction_type_last_interaction_class(
    simulation_rpacket_tracking,
):
    """Last interaction types of rpacket from LastInteractionTracker class"""
    interaction_type = (
        simulation_rpacket_tracking.transport.transport_state.last_interaction_type
    )
    return interaction_type


@pytest.fixture
def shell_id_last_interaction_class(
    simulation_rpacket_tracking,
):
    """
    shell_id when last interaction is line from LastInteractionTracker class
    """
    interaction_type = (
        simulation_rpacket_tracking.transport.transport_state.last_interaction_type
    )
    mask = interaction_type == InteractionType.LINE
    shell_id = (
        simulation_rpacket_tracking.transport.transport_state.last_line_interaction_shell_id
    )
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
def rpacket_tracker(simulation_rpacket_tracking):
    "RPacketTracker object from the simulation" ""
    rpacket_tracker = (
        simulation_rpacket_tracking.transport.transport_state.rpacket_tracker
    )
    return rpacket_tracker


@pytest.fixture(scope="module")
def last_interaction_type_rpacket_tracker(rpacket_tracker):
    """Last interaction types of rpacket from RPacketTracker class"""
    no_of_packets = len(rpacket_tracker)
    interaction_type = np.empty(no_of_packets, dtype=np.int64)

    for i in range(no_of_packets):
        # the last interaction is the second last element since the last element
        # correspond to reabsorbed/emission of the packet
        interaction_type[i] = rpacket_tracker[i].interaction_type[-2]

    return interaction_type


@pytest.fixture
def shell_id_rpacket_tracker(
    rpacket_tracker, last_interaction_type_rpacket_tracker
):
    """
    shell_id when last interaction is line from RPacketTracker class
    """
    no_of_packets = len(rpacket_tracker)
    shell_id = np.empty(no_of_packets, dtype=np.int64)

    for i in range(no_of_packets):
        shell_id[i] = rpacket_tracker[i].shell_id[-2]

    mask = last_interaction_type_rpacket_tracker == InteractionType.LINE
    last_line_interaction_shell_id = shell_id[mask]

    return last_line_interaction_shell_id


@pytest.fixture
def nu_rpacket_tracker(rpacket_tracker):
    """Output nu of rpacket from RPacketTracker class"""
    no_of_packets = len(rpacket_tracker)
    nu = np.empty(no_of_packets, dtype=np.float64)

    for i in range(no_of_packets):
        nu[i] = rpacket_tracker[i].nu[-2]

    return nu


def test_extend_array():
    rpacket_tracker = RPacketTracker(10)
    array = np.array([1, 2, 3, 4, 5], dtype=np.int64)

    new_array = rpacket_tracker.extend_array(array, array.size)

    assert new_array.size == array.size * rpacket_tracker.extend_factor
    assert new_array.dtype == array.dtype
    npt.assert_allclose(array, new_array[: array.size])


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
def test_rpacket_tracker_properties(expected, obtained, request):
    expected = request.getfixturevalue(expected)
    obtained = request.getfixturevalue(obtained)
    npt.assert_allclose(expected, obtained)


def test_boundary_interactions(rpacket_tracker, regression_data):
    no_of_packets = len(rpacket_tracker)

    max_boundary_interaction_size = max([tracker.boundary_interaction.size for tracker in rpacket_tracker])
    obtained_boundary_interaction = np.full(
        (no_of_packets, max_boundary_interaction_size),
        [-1],
        dtype=rpacket_tracker[0].boundary_interaction.dtype,
    )

    for i, tracker in enumerate(rpacket_tracker):
        obtained_boundary_interaction[
            i, : tracker.boundary_interaction.size
        ] = tracker.boundary_interaction

    expected_boundary_interaction = regression_data.sync_ndarray(
        obtained_boundary_interaction
    )
    npt.assert_array_equal(
        obtained_boundary_interaction, expected_boundary_interaction
    )


def test_rpacket_trackers_to_dataframe(simulation_rpacket_tracking):
    transport_state = simulation_rpacket_tracking.transport.transport_state
    rtracker_df = rpacket_trackers_to_dataframe(transport_state.rpacket_tracker)

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
