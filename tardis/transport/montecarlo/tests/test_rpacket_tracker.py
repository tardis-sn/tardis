import numpy.testing as npt
import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.base import run_tardis
from tardis.transport.montecarlo.packet_trackers import (
    rpacket_trackers_to_dataframe,
)


@pytest.fixture(scope="module")
def config_rpacket_tracker(example_configuration_dir):
    """Config object for rpacket tracker"""
    return Configuration.from_yaml(
        example_configuration_dir / "tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="module")
def simulation_rpacket_tracking_enabled(config_rpacket_tracker, atomic_dataset):
    """Simulation object with track_rpacket enabled"""
    config_rpacket_tracker.montecarlo.iterations = 3
    config_rpacket_tracker.montecarlo.no_of_packets = 4000
    config_rpacket_tracker.montecarlo.last_no_of_packets = -1
    config_rpacket_tracker.montecarlo.tracking.track_rpacket = True
    config_rpacket_tracker.spectrum.num = 2000
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config_verysimple,
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
    interation_type = simulation_rpacket_tracking_enabled.transport.transport_state.last_interaction_type
    return interaction_type


@pytest.fixture()
def shell_id_last_interaction_class(
    simulation_rpacket_tracking_enabled,
):
    """Last interaction types of rpacket from LastInteractionTracker"""
    interation_type = simulation_rpacket_tracking_enabled.transport.transport_state.last_interaction_type
    mask = interaction_type == 2
    shell_id = simulation_rpacket_tracking_enabled.transport.transport_state.last_line_interaction_shell_id
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


@pytest.fixture(scope="moduel")
def rpacket_tracker(simulation_rpacket_tracking_enabled):
    rpacket_tracker = simulation_rpacket_tracking_enabled.transport.transport_state.rpacket_tracker
    return rpacket_tracker


@pytest.fixture(scope="module")
def last_interaction_type_rpacket_tracker(rpacket_tracker):
    no_of_packets = len(rpacket_tracker)
    interaction_type = np.empty(no_of_packets, dtype = np.int64)

    for i in range(no_of_packets):
        interactions = rpacket_tracker[i].num_interactions
        interaction_type[i] = rpacket_tracker[i].interaction_type[interactions-1]

    return interaction_type


@pytest.fixture()
def shell_id_rpacket_tracker(rpacket_tracker, last_interaction_type_rpacket_tracker):
    no_of_packets = len(rpacket_tracker)
    shell_id = np.empty(no_of_packets, dtype = np.int64)

    for i in range(no_of_packets):
        interactions = rpacket_tracker[i].num_interactions
        shell_id[i] = rpacket_tracker[i].shell_id[interactions-1]
    mask = last_interaction_type_rpacket_tracker == 2
    last_line_interaction_shell_id = shell_id[mask]

    return last_line_interaction_shell_id


@pytest.fixture()
def nu_rpacket_tracker(rpacket_tracker):
    no_of_packets = len(rpacket_tracker)
    nu = np.empty(no_of_packets, dtype = np.float64)

    for i in range(no_of_packets):
        interactions = rpacket_tracker[i].num_interactions
        nu[i] = rpacket_tracker[i].nu[interactions-1]

    return nu
