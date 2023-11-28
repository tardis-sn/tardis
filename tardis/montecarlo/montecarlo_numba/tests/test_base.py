import pytest
import pandas as pd
import os
import numpy.testing as npt
from copy import deepcopy
from tardis.base import run_tardis
from pandas.testing import assert_frame_equal

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.montecarlo.montecarlo_numba.base import montecarlo_radial1d
from tardis.simulation import Simulation


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d(
    config_montecarlo_1e5_verysimple,
    atomic_dataset,
):
    atomic_data = deepcopy(atomic_dataset)
    config_montecarlo_1e5_verysimple.montecarlo.last_no_of_packets = 1e5
    config_montecarlo_1e5_verysimple.montecarlo.no_of_virtual_packets = 0
    config_montecarlo_1e5_verysimple.montecarlo.iterations = 1
    config_montecarlo_1e5_verysimple.plasma.line_interaction_type = "macroatom"

    montecarlo_configuration.ENABLE_VPACKET_TRACKING = True

    del config_montecarlo_1e5_verysimple["config_dirname"]

    sim = Simulation.from_config(
        config_montecarlo_1e5_verysimple, atom_data=atomic_data
    )
    transport_state = sim.transport.initialize_transport_state(
        sim.simulation_state,
        sim.plasma,
        no_of_packets=int(1e5),
        no_of_virtual_packets=3,
        iteration=0,
    )
    montecarlo_radial1d(
        transport_state,
        sim.simulation_state.time_explosion,
        iteration=0,
        total_iterations=5,
        show_progress_bars=False,
    )

    assert transport_state.no_of_packets == 1e5


def test_montecarlo_main_loop(
    config_montecarlo_1e5_verysimple,
    atomic_dataset,
    tardis_ref_path,
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
    sim.run_convergence()
    sim.run_final()

    compare_fname = os.path.join(
        tardis_ref_path, "montecarlo_1e5_compare_data.h5"
    )
    if request.config.getoption("--generate-reference"):
        sim.to_hdf(compare_fname, overwrite=True)

    # Load compare data from refdata
    expected_nu = pd.read_hdf(
        compare_fname, key="/simulation/transport/output_nu"
    ).values
    expected_energy = pd.read_hdf(
        compare_fname, key="/simulation/transport/output_energy"
    ).values
    expected_nu_bar_estimator = pd.read_hdf(
        compare_fname, key="/simulation/transport/nu_bar_estimator"
    ).values
    expected_j_estimator = pd.read_hdf(
        compare_fname, key="/simulation/transport/j_estimator"
    ).values

    actual_energy = (
        sim.transport.transport_state.packet_collection.output_energies
    )
    actual_nu = sim.transport.transport_state.packet_collection.output_nus
    actual_nu_bar_estimator = (
        sim.transport.transport_state.estimators.nu_bar_estimator
    )
    actual_j_estimator = sim.transport.transport_state.estimators.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-13)
