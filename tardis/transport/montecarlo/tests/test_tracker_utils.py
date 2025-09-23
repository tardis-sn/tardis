import numpy as np
from numba import typeof

from tardis.transport.montecarlo.packets.RPacketLastInteractionTracker import RPacketLastInteractionTracker
from tardis.transport.montecarlo.packets.trackers.tracker_full import TrackerFull
from tardis.transport.montecarlo.packets.packet_trackers import (
    generate_rpacket_last_interaction_tracker_list,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full_util import generate_tracker_full_list


def test_generate_rpacket_tracker_list():
    no_of_packets = 10
    length = 10
    random_index = np.random.randint(0, no_of_packets)

    rpacket_tracker_list = generate_tracker_full_list(no_of_packets, length)

    assert len(rpacket_tracker_list) == no_of_packets
    assert len(rpacket_tracker_list[random_index].shell_id) == length
    assert typeof(rpacket_tracker_list[random_index]) == typeof(
        TrackerFull(length)
    )


def test_generate_rpacket_last_interaction_tracker_list():
    no_of_packets = 50
    random_index = np.random.randint(0, no_of_packets)

    rpacket_last_interaction_tracker_list = (
        generate_rpacket_last_interaction_tracker_list(no_of_packets)
    )

    assert len(rpacket_last_interaction_tracker_list) == no_of_packets
    assert typeof(
        rpacket_last_interaction_tracker_list[random_index]
    ) == typeof(RPacketLastInteractionTracker())
