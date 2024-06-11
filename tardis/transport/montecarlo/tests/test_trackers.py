import pytest
import numpy as np
import numpy.testing as npt

from tardis.transport.montecarlo.packet_trackers import (
    RPacketLastInteractionTracker,
)


@pytest.fixture()
def last_interactions(simulation_rpacket_tracking_enabled):
    sim = simulation_rpacket_tracking_enabled
    transport_state = sim.transport.transport_state
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


def test_tracking_from_last_interaction_simulation(last_interactions):
    checkedTypes, toCheckTypes = last_interactions
    npt.assert_array_equal(checkedTypes, toCheckTypes)


def test_tracking_from_rpacket_tracker_simulation(
    simulation_rpacket_tracking_enabled,
):
    sim = simulation_rpacket_tracking_enabled
    rpacket_trackers = sim.transport.transport_state.rpacket_tracker
    rpacket_last_interaction_trackers = (
        sim.transport.transport_state.rpacket_last_interaction_trackers
    )
    no_of_packets = len(rpacket_trackers)

    index_from_rpacket_tracker = np.empty(no_of_packets, dtype=np.int64)
    index_from_last_interaction_tracker = np.empty(
        no_of_packets, dtype=np.int64
    )

    r_from_rpacket_tracker = np.empty(no_of_packets, dtype=np.float64)
    r_from_last_interaction_tracker = np.empty(no_of_packets, dtype=np.float64)

    nu_from_rpacket_tracker = np.empty(no_of_packets, dtype=np.float64)
    nu_from_last_interaction_tracker = np.empty(no_of_packets, dtype=np.float64)

    energy_from_rpacket_tracker = np.empty(no_of_packets, dtype=np.float64)
    energy_from_last_interaction_tracker = np.empty(
        no_of_packets, dtype=np.float64
    )

    shell_id_from_rpacket_tracker = np.empty(no_of_packets, dtype=np.int64)
    shell_id_from_last_interaction_tracker = np.empty(
        no_of_packets, dtype=np.int64
    )

    type_from_rpacket_tracker = np.empty(no_of_packets, dtype=np.int64)
    type_from_last_interaction_tracker = np.empty(no_of_packets, dtype=np.int64)

    for i in range(no_of_packets):
        interactions = rpacket_trackers[i].num_interactions

        index_from_rpacket_tracker[i] = rpacket_trackers[i].index
        index_from_last_interaction_tracker[
            i
        ] = rpacket_last_interaction_trackers[i].index

        r_from_rpacket_tracker[i] = rpacket_trackers[i].r[interactions - 1]
        r_from_last_interaction_tracker[i] = rpacket_last_interaction_trackers[
            i
        ].r

        nu_from_rpacket_tracker[i] = rpacket_trackers[i].nu[interactions - 1]
        nu_from_last_interaction_tracker[i] = rpacket_last_interaction_trackers[
            i
        ].nu

        energy_from_rpacket_tracker[i] = rpacket_trackers[i].energy[
            interactions - 1
        ]
        energy_from_last_interaction_tracker[
            i
        ] = rpacket_last_interaction_trackers[i].energy

        shell_id_from_rpacket_tracker[i] = rpacket_trackers[i].shell_id[
            interactions - 1
        ]
        shell_id_from_last_interaction_tracker[
            i
        ] = rpacket_last_interaction_trackers[i].shell_id

        type_from_rpacket_tracker[i] = rpacket_trackers[i].interaction_type[
            interactions - 1
        ]
        type_from_last_interaction_tracker[
            i
        ] = rpacket_last_interaction_trackers[i].interaction_type

    npt.assert_array_equal(
        index_from_rpacket_tracker, index_from_last_interaction_tracker
    )
    npt.assert_allclose(r_from_rpacket_tracker, r_from_last_interaction_tracker)
    npt.assert_allclose(
        nu_from_rpacket_tracker, nu_from_last_interaction_tracker
    )
    npt.assert_allclose(
        energy_from_rpacket_tracker,
        energy_from_last_interaction_tracker,
    )
    npt.assert_array_equal(
        shell_id_from_rpacket_tracker, shell_id_from_last_interaction_tracker
    )
    npt.assert_array_equal(
        type_from_rpacket_tracker, type_from_last_interaction_tracker
    )
