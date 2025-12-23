import numpy as np
import numpy.testing as npt
import pytest

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType
from tardis.transport.montecarlo.packets.trackers.tracker_last_interaction import TrackerLastInteraction

@pytest.fixture(scope="module")
def interaction_type_in_use(
    nb_simulation_verysimple,
):
    """Last interaction types of rpacket from LastInteractionTracker class"""
    transport_state = nb_simulation_verysimple.transport.transport_state
    df = transport_state.tracker_last_interaction_df
    return df['last_interaction_type'].astype(str).values


@pytest.fixture
def shell_id_in_use(
    nb_simulation_verysimple,
    interaction_type_in_use,
):
    """
    `shell_id` when last interaction is line from LastInteractionTracker class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    df = transport_state.tracker_last_interaction_df
    shell_id = df['shell_id'].values
    mask = interaction_type_in_use == "LINE"
    return shell_id[mask]


@pytest.fixture
def r_in_use(
    nb_simulation_verysimple,
    interaction_type_in_use,
):
    """
    `r` when last interaction is line from LastInteractionTracker class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    df = transport_state.tracker_last_interaction_df
    r = df['radius'].values
    mask = interaction_type_in_use == "LINE"
    return r[mask]


@pytest.fixture(scope="module")
def interaction_type_to_check(
    nb_simulation_verysimple,
):
    """
    Last interaction types of rpacket from TrackerLastInteraction class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    df = transport_state.tracker_last_interaction_df
    return df['last_interaction_type'].astype(str).values


@pytest.fixture
def shell_id_to_check(
    nb_simulation_verysimple,
    interaction_type_to_check,
):
    """
    shell_id when last interaction is line from TrackerLastInteraction class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    df = transport_state.tracker_last_interaction_df
    shell_id = df['shell_id'].values
    mask = interaction_type_to_check == "LINE"
    return shell_id[mask]


@pytest.fixture
def r_to_check(
    nb_simulation_verysimple,
    interaction_type_to_check,
):
    """
    `r` when last interaction is line from TrackerLastInteraction class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    df = transport_state.tracker_last_interaction_df
    r = df['radius'].values
    mask = interaction_type_to_check == "LINE"
    return r[mask]


@pytest.fixture
def nu_packet_collection(
    nb_simulation_verysimple,
):
    """Last interaction output nus of rpacket from packet_collection"""
    transport_state = nb_simulation_verysimple.transport.transport_state
    df = transport_state.tracker_last_interaction_df
    return df['after_nu'].values


@pytest.fixture
def nu_to_check(
    nb_simulation_verysimple,
):
    """
    Last interaction output nus of rpacket from TrackerLastInteraction class
    """
    transport_state = nb_simulation_verysimple.transport.transport_state
    df = transport_state.tracker_last_interaction_df
    return df['after_nu'].values


def test_defaults():
    tracker = TrackerLastInteraction()
    assert tracker.shell_id == -1
    assert tracker.interaction_type == -1
    assert tracker.interactions_count == 0
    npt.assert_almost_equal(tracker.r, -1.0)
    assert tracker.interaction_line_absorb_id == -1
    assert tracker.interaction_line_emit_id == -1


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
    npt.assert_array_equal(expected, obtained)
