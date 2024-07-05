from copy import deepcopy

import numpy as np
import numpy.testing as npt
import pytest

from tardis.base import run_tardis
from tardis.transport.montecarlo.packet_trackers import (
    rpacket_trackers_to_dataframe,
)


@pytest.fixture(scope="module")
def simulation_rpacket_tracking_enabled(config_verysimple, atomic_dataset):
    config_verysimple.montecarlo.iterations = 3
    config_verysimple.montecarlo.no_of_packets = 4000
    config_verysimple.montecarlo.last_no_of_packets = -1
    config_verysimple.spectrum.virtual.virtual_packet_logging = True
    config_verysimple.montecarlo.no_of_virtual_packets = 1
    config_verysimple.montecarlo.tracking.track_rpacket = True
    config_verysimple.spectrum.num = 2000
    atomic_data = deepcopy(atomic_dataset)
    sim = run_tardis(
        config_verysimple,
        atom_data=atomic_data,
        show_convergence_plots=False,
    )
    return sim


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
