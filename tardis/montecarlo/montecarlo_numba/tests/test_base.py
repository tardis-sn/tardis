import pytest
import pandas as pd
import os
import numpy.testing as npt
import numpy as np
from copy import deepcopy
from tardis.base import run_tardis
from pandas.testing import assert_frame_equal

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.simulation import Simulation
from tardis.montecarlo.montecarlo_numba.base import (
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


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d():
    assert False


def test_rpacket_trackers_to_dataframe(simulation_rpacket_tracking_enabled):
    sim = simulation_rpacket_tracking_enabled
    rtracker_df = rpacket_trackers_to_dataframe(sim.runner.rpacket_tracker)

    # check df shape and column names
    assert rtracker_df.shape == (
        sum([len(tracker.r) for tracker in sim.runner.rpacket_tracker]),
        6,
    )
    assert npt.assert_array_equal(
        sim.runner.rpacket_tracker_df.columns.values,
        ["status", "r", "nu", "mu", "energy", "shell_id"],
    )

    # check all data with rpacket_tracker
    expected_rtrackers = []
    for rpacket in sim.runner.rpacket_tracker:
        for rpacket_step_no in range(len(rpacket.r)):
            expected_rtrackers.append(
                [
                    rpacket.status[rpacket_step_no],
                    rpacket.r[rpacket_step_no],
                    rpacket.nu[rpacket_step_no],
                    rpacket.mu[rpacket_step_no],
                    rpacket.energy[rpacket_step_no],
                    rpacket.shell_id[rpacket_step_no],
                ]
            )
    npt.assert_array_equal(rtracker_df.to_numpy(), np.array(expected_rtrackers))


def test_montecarlo_main_loop(
    config_montecarlo_1e5_verysimple,
    atomic_dataset,
    tardis_ref_path,
    tmpdir,
    set_seed_fixture,
    random_call_fixture,
    request,
):

    montecarlo_configuration.LEGACY_MODE_ENABLED = True
    # Setup model config from verysimple
    atomic_data = deepcopy(atomic_dataset)
    config_montecarlo_1e5_verysimple.montecarlo.last_no_of_packets = 1e5
    config_montecarlo_1e5_verysimple.montecarlo.no_of_virtual_packets = 0
    config_montecarlo_1e5_verysimple.montecarlo.iterations = 1
    config_montecarlo_1e5_verysimple.plasma.line_interaction_type = "macroatom"
    del config_montecarlo_1e5_verysimple["config_dirname"]

    sim = Simulation.from_config(
        config_montecarlo_1e5_verysimple, atom_data=atomic_data
    )
    sim.run()

    compare_fname = os.path.join(
        tardis_ref_path, "montecarlo_1e5_compare_data.h5"
    )
    if request.config.getoption("--generate-reference"):
        sim.to_hdf(compare_fname, overwrite=True)

    # Load compare data from refdata
    expected_nu = pd.read_hdf(
        compare_fname, key="/simulation/runner/output_nu"
    ).values
    expected_energy = pd.read_hdf(
        compare_fname, key="/simulation/runner/output_energy"
    ).values
    expected_nu_bar_estimator = pd.read_hdf(
        compare_fname, key="/simulation/runner/nu_bar_estimator"
    ).values
    expected_j_estimator = pd.read_hdf(
        compare_fname, key="/simulation/runner/j_estimator"
    ).values

    actual_energy = sim.runner.output_energy
    actual_nu = sim.runner.output_nu
    actual_nu_bar_estimator = sim.runner.nu_bar_estimator
    actual_j_estimator = sim.runner.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy.value, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu.value, expected_nu, rtol=1e-13)
