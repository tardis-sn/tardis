import pytest
import numpy as np
import numpy.testing as npt

from tardis.transport.montecarlo.packet_trackers import (
    RPacketLastInteractionTracker,
)


@pytest.fixture()
def interaction_type_last_interaction_class1(
    nb_simulation_verysimple,
):
    """Last interaction types of rpacket from LastInteractionTracker"""
    sim = nb_simulation_verysimple
    transport_state = sim.transport.transport_state
    interaction_type = transport_state.last_interaction_type
    return interaction_type


@pytest.fixture()
def interaction_type_last_interaction_class2(
    nb_simulation_verysimple,
):
    """Last interaction types of rpacket from RPacketLastInteractionTracker"""
    sim = nb_simulation_verysimple
    transport_state = sim.transport.transport_state
    interaction_type = np.empty(
        len(transport_state.rpacket_tracker), dtype=np.int64
    )
    for i, last_interaction_tracker in enumerate(
        transport_state.rpacket_tracker
    ):
        interaction_type[i] = last_interaction_tracker.interaction_type

    return interaction_type


@pytest.fixture()
def nu_from_packet_collection(
    nb_simulation_verysimple,
):
    """Last interaction output nus of rpacket from packet_collection"""
    sim = nb_simulation_verysimple
    packet_collection = sim.transport.transport_state.packet_collection
    return packet_collection.output_nus


@pytest.fixture()
def nu_from_last_interaction_class(
    nb_simulation_verysimple,
):
    """Last interaction output nus of rpacket from RPacketLastInteractionTracker"""
    sim = nb_simulation_verysimple
    transport_state = sim.transport.transport_state
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


def test_last_interaction_type(
    interaction_type_last_interaction_class1,
    interaction_type_last_interaction_class2,
):
    npt.assert_array_equal(
        interaction_type_last_interaction_class1,
        interaction_type_last_interaction_class2,
    )


def test_last_interaction_nu(
    nu_from_packet_collection,
    nu_from_last_interaction_class,
):
    npt.assert_allclose(
        nu_from_packet_collection, nu_from_last_interaction_class
    )
