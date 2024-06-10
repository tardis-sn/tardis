from tardis.transport.montecarlo.packet_trackers import (
    RPacketLastInteractionTracker,
)
from tardis.transport.montecarlo.r_packet import RPacket


def test_defaults():
    tracker = RPacketLastInteractionTracker()
    assert tracker.index == -1
    assert tracker.r == -1
    assert tracker.nu == 0
    assert tracker.energy == 0
    assert tracker.shell_id == -1
    assert tracker.interaction_type == -1
