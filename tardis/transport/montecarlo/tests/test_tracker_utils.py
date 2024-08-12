import pytest
import numpy as np
from numba import typeof

from tardis.transport.montecarlo.packet_trackers import (
    RPacketTracker,
    RPacketLastInteractionTracker,
    generate_rpacket_tracker_list,
    generate_rpacket_last_interaction_tracker_list,
)


def test_generate_rpacket_tracker_list():
    no_of_packets = 10
    length = 10
    random_index = np.random.randint(0, no_of_packets)

    rpacket_tracker_list = generate_rpacket_tracker_list(no_of_packets, length)

    assert len(rpacket_tracker_list) == no_of_packets
    assert len(rpacket_tracker_list[random_index].shell_id) == length
    assert typeof(rpacket_tracker_list[random_index]) == typeof(
        RPacketTracker(length)
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
