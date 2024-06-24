import pytest
import numpy as np
import numpy.testing as npt

from tardis.transport.montecarlo.packet_trackers import (
    RPacketLastInteractionTracker,
)


@pytest.fixture(scope="module")
def interaction_type_last_interaction_class_old(
    nb_simulation_verysimple,
):
    """Last interaction types of rpacket from LastInteractionTracker"""
    transport_state = nb_simulation_verysimple.transport.transport_state
    interaction_type = transport_state.last_interaction_type
    return interaction_type


@pytest.fixture(scope="module")
def interaction_type_last_interaction_class_new(
    nb_simulation_verysimple,
):
    """Last interaction types of rpacket from RPacketLastInteractionTracker"""
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
def last_line_interaction_old(
    nb_simulation_verysimple,
    interaction_type_last_interaction_class_old,
):
    """Last line interaction data of rpacket from LastInteractionTracker"""
    transport_state = nb_simulation_verysimple.transport.transport_state
    mask = interaction_type_last_interaction_class_old == 2

    def get_attribute_data(attribute):
        attribute_data = getattr(transport_state, attribute)
        attribute_data = attribute_data[mask]
        return attribute_data

    return get_attribute_data


@pytest.fixture()
def last_line_interaction_new(
    nb_simulation_verysimple,
    interaction_type_last_interaction_class_new,
):
    """Last line interaction data of rpacket from RPacketLastInteractionTracker"""
    rpacket_tracker = (
        nb_simulation_verysimple.transport.transport_state.rpacket_tracker
    )
    mask = interaction_type_last_interaction_class_new == 2

    def get_attribute_data(attribute):
        attribute_data = np.empty(len(rpacket_tracker))
        for i, last_interaction_tracker in enumerate(rpacket_tracker):
            attribute_data[i] = getattr(last_interaction_tracker, attribute)
        attribute_data = attribute_data[mask]
        return attribute_data

    return get_attribute_data


@pytest.fixture()
def nu_from_packet_collection(
    nb_simulation_verysimple,
):
    """Last interaction output nus of rpacket from packet_collection"""
    packet_collection = (
        nb_simulation_verysimple.transport.transport_state.packet_collection
    )
    return packet_collection.output_nus


@pytest.fixture()
def nu_from_last_interaction_class_new(
    nb_simulation_verysimple,
):
    """Last interaction output nus of rpacket from RPacketLastInteractionTracker"""
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
            "interaction_type_last_interaction_class_old",
            "interaction_type_last_interaction_class_new",
        ),
        ("nu_from_packet_collection", "nu_from_last_interaction_class_new"),
    ],
)
def test_last_interaction_properties(expected, obtained, request):
    expected = request.getfixturevalue(expected)
    obtained = request.getfixturevalue(obtained)
    npt.assert_allclose(expected, obtained)


@pytest.mark.parametrize(
    "argument_old,argument_new",
    [
        (
            "last_interaction_in_nu",
            "interaction_in_line_nu",
        ),
        (
            "last_line_interaction_in_id",
            "interaction_in_line_id",
        ),
        (
            "last_line_interaction_out_id",
            "interaction_out_line_id",
        ),
        (
            "last_line_interaction_shell_id",
            "shell_id",
        ),
    ],
)
def test_last_line_interaction_properties(
    last_line_interaction_old,
    last_line_interaction_new,
    argument_old,
    argument_new,
):
    npt.assert_allclose(
        last_line_interaction_old(argument_old),
        last_line_interaction_new(argument_new),
    )
