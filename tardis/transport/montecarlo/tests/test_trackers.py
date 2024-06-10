import pytest
import numpy as np

from tardis.transport.montecarlo.packet_trackers import (
    RPacketLastInteractionTracker,
)


@pytest.fixture()
def last_interactions(nb_simulation_verysimple):
    transport_state = nb_simulation_verysimple.transport.transport_state
    checkedTypes = transport_state.last_interaction_type
    toCheckTypes = np.empty(
        len(transport_state.rpacket_last_interaction_trackers), dtype=np.int64
    )
    for i, rpacket_last_interaction_tracker in enumerate(
        transport_state.rpacket_last_interaction_trackers
    ):
        toCheckTypes[i] = rpacket_last_interaction_tracker.interaction_type

    return (checkedTypes, toCheckTypes)


def test_defaults():
    tracker = RPacketLastInteractionTracker()
    assert tracker.index == -1
    assert tracker.r == -1
    assert tracker.nu == 0
    assert tracker.energy == 0
    assert tracker.shell_id == -1
    assert tracker.interaction_type == -1


def test_tracking_manual(static_packet):
    tracker = RPacketLastInteractionTracker()
    tracker.track_last_interaction(static_packet)
    assert tracker.index == 0
    assert tracker.r == 7.5e14
    assert tracker.nu == 0.4
    assert tracker.energy == 0.9


def test_tracking_simulation(last_interactions):
    checkedTypes, toCheckTypes = last_interactions
    assert np.array_equal(checkedTypes, toCheckTypes)
