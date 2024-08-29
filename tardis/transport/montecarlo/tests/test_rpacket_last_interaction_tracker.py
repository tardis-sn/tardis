import pytest
import numpy as np
import numpy.testing as npt

from tardis.transport.montecarlo.r_packet import InteractionType
from tardis.transport.montecarlo.packet_trackers import (
    RPacketLastInteractionTracker,
)


@pytest.fixture(scope="module")
def interaction_type_in_use(
    nb_simulation_verysimple,
):
    """Last interaction types of rpacket from LastInteractionTracker class"""
    transport_state = nb_simulation_verysimple.transport.transport_state
    interaction_type = transport_state.last_interaction_type
    return interaction_type


@pytest.fixture()
def shell_id_in_use(
    nb_simulation_verysimple,
    interaction_type_in_use,
):
    """
    `shell_id` when last interaction is line from LastInteractionTracker class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    shell_id = transport_state.last_line_interaction_shell_id
    mask = interaction_type_in_use == InteractionType.LINE
    return shell_id[mask]


@pytest.fixture()
def r_in_use(
    nb_simulation_verysimple,
    interaction_type_in_use,
):
    """
    `r` when last interaction is line from LastInteractionTracker class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    r = transport_state.last_interaction_in_r
    mask = interaction_type_in_use == InteractionType.LINE
    return r[mask]


@pytest.fixture(scope="module")
def interaction_type_to_check(
    nb_simulation_verysimple,
):
    """
    Last interaction types of rpacket from RPacketLastInteractionTracker class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    interaction_type = np.empty(
        len(transport_state.rpacket_tracker), dtype=np.int64
    )
    for i, last_interaction_tracker in enumerate(
        transport_state.rpacket_tracker
    ):
        interaction_type[i] = last_interaction_tracker.interaction_type

    return interaction_type


@pytest.fixture()
def shell_id_to_check(
    nb_simulation_verysimple,
    interaction_type_to_check,
):
    """
    shell_id when last interaction is line from RPacketLastInteractionTracker class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    shell_id = np.empty(len(transport_state.rpacket_tracker), dtype=np.int64)
    for i, last_interaction_tracker in enumerate(
        transport_state.rpacket_tracker
    ):
        shell_id[i] = last_interaction_tracker.shell_id
    mask = interaction_type_to_check == InteractionType.LINE
    return shell_id[mask]


@pytest.fixture()
def r_to_check(
    nb_simulation_verysimple,
    interaction_type_to_check,
):
    """
    `r` when last interaction is line from RPacketLastInteractionTracker class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    r = np.empty(len(transport_state.rpacket_tracker), dtype=np.int64)
    for i, last_interaction_tracker in enumerate(
        transport_state.rpacket_tracker
    ):
        r[i] = last_interaction_tracker.r
    mask = interaction_type_to_check == InteractionType.LINE
    return r[mask]


@pytest.fixture()
def nu_packet_collection(
    nb_simulation_verysimple,
):
    """Last interaction output nus of rpacket from packet_collection"""
    packet_collection = (
        nb_simulation_verysimple.transport.transport_state.packet_collection
    )
    return packet_collection.output_nus


@pytest.fixture()
def nu_to_check(
    nb_simulation_verysimple,
):
    """
    Last interaction output nus of rpacket from RPacketLastInteractionTracker class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    nu = np.empty(len(transport_state.rpacket_tracker), dtype=np.float64)
    for i, last_interaction_tracker in enumerate(
        transport_state.rpacket_tracker
    ):
        nu[i] = last_interaction_tracker.nu

    return nu


def test_defaults():
    tracker = RPacketLastInteractionTracker()
    assert tracker.index == -1
    assert tracker.shell_id == -1
    assert tracker.interaction_type == -1
    npt.assert_almost_equal(tracker.r, -1.0)
    npt.assert_almost_equal(tracker.nu, 0.0)
    npt.assert_almost_equal(tracker.energy, 0.0)


def test_tracking_manual(static_packet):
    tracker = RPacketLastInteractionTracker()
    tracker.track(static_packet)
    assert tracker.index == 0
    npt.assert_almost_equal(tracker.r, 7.5e14)
    npt.assert_almost_equal(tracker.nu, 0.4)
    npt.assert_almost_equal(tracker.energy, 0.9)


@pytest.mark.parametrize(
    "expected,obtained",
    [
        (
            "interaction_type_in_use",
            "interaction_type_to_check",
        ),
        (
            "shell_id_in_use",
            "shell_id_to_check",
        ),
        (
            "r_in_use",
            "r_to_check",
        ),
        ("nu_packet_collection", "nu_to_check"),
    ],
)
def test_last_interaction_properties(expected, obtained, request):
    expected = request.getfixturevalue(expected)
    obtained = request.getfixturevalue(obtained)
    npt.assert_allclose(expected, obtained)
