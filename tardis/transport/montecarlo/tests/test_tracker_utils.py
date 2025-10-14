import numpy as np
from numba import typeof

from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import TrackerLastInteraction
from tardis.transport.montecarlo.packets.trackers.tracker_full import TrackerFull
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction_util import (
    generate_tracker_last_interaction_list,
)
from tardis.transport.montecarlo.packets.trackers.tracker_full_util import generate_tracker_full_list


def test_generate_rpacket_tracker_list():
    no_of_packets = 10
    length = 10
    extend_factor = 2
    random_index = np.random.randint(0, no_of_packets)

    rpacket_tracker_list = generate_tracker_full_list(no_of_packets, length)

    assert len(rpacket_tracker_list) == no_of_packets
    assert len(rpacket_tracker_list[random_index].before_shell_id) == length
    assert typeof(rpacket_tracker_list[random_index]) == typeof(
        TrackerFull(length, extend_factor)
    )


def test_generate_rpacket_last_interaction_tracker_list():
    no_of_packets = 50
    random_index = np.random.randint(0, no_of_packets)

    rpacket_last_interaction_tracker_list = (
        generate_tracker_last_interaction_list(no_of_packets)
    )

    assert len(rpacket_last_interaction_tracker_list) == no_of_packets
    assert typeof(
        rpacket_last_interaction_tracker_list[random_index]
    ) == typeof(TrackerLastInteraction())
