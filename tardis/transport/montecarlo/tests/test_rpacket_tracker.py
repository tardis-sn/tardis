from copy import deepcopy

import numpy.testing as npt
import pytest
import numpy as np

from tardis.io.configuration.config_reader import Configuration
from tardis.transport.montecarlo.r_packet import InteractionType
from tardis.base import run_tardis


@pytest.fixture(scope="module")
def config_rpacket_tracker(example_configuration_dir):
    """Config object for rpacket tracker"""
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="module")
def simulation_rpacket_tracking_enabled(config_rpacket_tracker, atomic_dataset):
    """Simulation object with track_rpacket_enabled"""
    config_rpacket_tracker.montecarlo.iterations = 1
    config_rpacket_tracker.montecarlo.no_of_packets = 4000
    config_rpacket_tracker.montecarlo.last_no_of_packets = -1
    config_rpacket_tracker.montecarlo.tracking.track_rpacket = True
    config_rpacket_tracker.spectrum.num = 2000
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config_rpacket_tracker,
        atom_data=atomic_data,
        show_convergence_plots=False,
    )
    return sim


@pytest.fixture()
def interaction_type_last_interaction_class(
    simulation_rpacket_tracking_enabled,
):
    """Last interaction types of rpacket from LastInteractionTracker"""
    interaction_type = (
        simulation_rpacket_tracking_enabled.transport.transport_state.last_interaction_type
    )
    return interaction_type


@pytest.fixture()
def shell_id_last_interaction_class(
    simulation_rpacket_tracking_enabled,
):
    """Last Line Interaction Shell Id of rpacket from LastInteractionTracker"""
    interaction_type = (
        simulation_rpacket_tracking_enabled.transport.transport_state.last_interaction_type
    )
    mask = interaction_type == InteractionType.LINE
    shell_id = (
        simulation_rpacket_tracking_enabled.transport.transport_state.last_line_interaction_shell_id
    )
    last_line_interaction_shell_id = shell_id[mask]

    return last_line_interaction_shell_id


@pytest.fixture()
def nu_from_packet_collection(
    simulation_rpacket_tracking_enabled,
):
    """Last interaction output nus of rpacket from packet_collection"""
    packet_collection = (
        simulation_rpacket_tracking_enabled.transport.transport_state.packet_collection
    )
    return packet_collection.output_nus


@pytest.fixture(scope="module")
def rpacket_tracker(simulation_rpacket_tracking_enabled):
    "RPacketTracker object from the simulation" ""
    rpacket_tracker = (
        simulation_rpacket_tracking_enabled.transport.transport_state.rpacket_tracker
    )
    return rpacket_tracker


@pytest.fixture(scope="module")
def last_interaction_type_rpacket_tracker(rpacket_tracker):
    """Last interaction types of rpacket from RPacketTracker"""
    no_of_packets = len(rpacket_tracker)
    interaction_type = np.empty(no_of_packets, dtype=np.int64)

    for i in range(no_of_packets):
        interactions = rpacket_tracker[i].num_interactions
        interaction_type[i] = rpacket_tracker[i].interaction_type[
            interactions - 1
        ]

    return interaction_type


@pytest.fixture()
def shell_id_rpacket_tracker(
    rpacket_tracker, last_interaction_type_rpacket_tracker
):
    """Last line interaction shell id of rpacket from RPacketTracker"""
    no_of_packets = len(rpacket_tracker)
    shell_id = np.empty(no_of_packets, dtype=np.int64)

    for i in range(no_of_packets):
        interactions = rpacket_tracker[i].num_interactions
        shell_id[i] = rpacket_tracker[i].shell_id[interactions - 1]
    mask = last_interaction_type_rpacket_tracker == InteractionType.LINE
    last_line_interaction_shell_id = shell_id[mask]

    return last_line_interaction_shell_id


@pytest.fixture()
def nu_rpacket_tracker(rpacket_tracker):
    """Output nu of rpacket from RPacketTracker"""
    no_of_packets = len(rpacket_tracker)
    nu = np.empty(no_of_packets, dtype=np.float64)

    for i in range(no_of_packets):
        interactions = rpacket_tracker[i].num_interactions
        nu[i] = rpacket_tracker[i].nu[interactions - 1]

    return nu


@pytest.mark.parametrize(
    "expected,obtained",
    [
        (
            "interaction_type_last_interaction_class",
            "last_interaction_type_rpacket_tracker",
        ),
        ("shell_id_last_interaction_class", "shell_id_rpacket_tracker"),
        ("nu_from_packet_collection", "nu_rpacket_tracker"),
    ],
)
def test_rpacket_tracker_properties(expected, obtained, request):
    expected = request.getfixturevalue(expected)
    obtained = request.getfixturevalue(obtained)
    npt.assert_allclose(expected, obtained)


def test_rpacket_trackers_to_dataframe(simulation_rpacket_tracking_enabled):
    sim = simulation_rpacket_tracking_enabled
    transport_state = sim.transport.transport_state
    rtracker_df = rpacket_trackers_to_dataframe(transport_state.rpacket_tracker)

    # check df shape and column names
    assert rtracker_df.shape == (
        sum([len(tracker.r) for tracker in transport_state.rpacket_tracker]),
        8,
    )
    npt.assert_array_equal(
        transport_state.rpacket_tracker_df.columns.values,
        np.array(
            [
                "status",
                "seed",
                "r",
                "nu",
                "mu",
                "energy",
                "shell_id",
                "interaction_type",
            ]
        ),
    )

    # check all data with rpacket_tracker
    expected_rtrackers = []
    for rpacket in transport_state.rpacket_tracker:
        for rpacket_step_no in range(len(rpacket.r)):
            expected_rtrackers.append(
                [
                    rpacket.status[rpacket_step_no],
                    rpacket.seed,
                    rpacket.r[rpacket_step_no],
                    rpacket.nu[rpacket_step_no],
                    rpacket.mu[rpacket_step_no],
                    rpacket.energy[rpacket_step_no],
                    rpacket.shell_id[rpacket_step_no],
                    rpacket.interaction_type[rpacket_step_no],
                ]
            )
    npt.assert_array_equal(rtracker_df.to_numpy(), np.array(expected_rtrackers))
