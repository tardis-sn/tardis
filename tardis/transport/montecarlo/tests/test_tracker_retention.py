from copy import deepcopy

import pytest

from tardis.io.atom_data import AtomData
from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation

NUMBER_OF_RPACKETS = 10


@pytest.fixture
def tracker_simulation(
    config_verysimple: Configuration, atomic_dataset: AtomData
) -> Simulation:
    """Very simple simulation that requests virtual packets but no tracking."""
    config = deepcopy(config_verysimple)
    config.montecarlo.no_of_virtual_packets = 3
    config.montecarlo.tracking.track_rpacket = False
    return Simulation.from_config(config, atom_data=deepcopy(atomic_dataset))


def test_record_full_histories_only_for_requested_run(
    tracker_simulation: Simulation,
) -> None:
    transport = tracker_simulation.transport

    tracker_simulation.iterate(NUMBER_OF_RPACKETS)
    assert transport.transport_state.tracker_full_df is None

    tracker_simulation.iterate(NUMBER_OF_RPACKETS, record_full_histories=True)
    tracker_full_df = transport.transport_state.tracker_full_df
    assert tracker_full_df.index.get_level_values("packet_id").nunique() == (
        NUMBER_OF_RPACKETS
    )
    assert (
        len(transport.transport_state.tracker_last_interaction_df)
        == NUMBER_OF_RPACKETS
    )

    tracker_simulation.iterate(NUMBER_OF_RPACKETS)
    assert transport.transport_state.tracker_full_df is None
