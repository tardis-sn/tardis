import pytest
import pandas as pd

import numpy.testing as npt
from copy import deepcopy
from tardis.base import run_tardis
from pandas.testing import assert_frame_equal

from tardis.montecarlo import (
    montecarlo_configuration as montecarlo_configuration,
)
from tardis.simulation import Simulation


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d():
    assert False


@pytest.fixture(scope="function")
def montecarlo_main_loop_config(
    config_montecarlo_1e5_verysimple,
):
    montecarlo_configuration.LEGACY_MODE_ENABLED = True
    # Setup model config from verysimple

    config_montecarlo_1e5_verysimple.montecarlo.last_no_of_packets = 1e5
    config_montecarlo_1e5_verysimple.montecarlo.no_of_virtual_packets = 0
    config_montecarlo_1e5_verysimple.montecarlo.iterations = 1
    config_montecarlo_1e5_verysimple.plasma.line_interaction_type = "macroatom"

    del config_montecarlo_1e5_verysimple["config_dirname"]
    return config_montecarlo_1e5_verysimple


def test_montecarlo_main_loop(
    montecarlo_main_loop_config,
    regression_data,
    atomic_dataset,
):
    atomic_dataset = deepcopy(atomic_dataset)
    montecarlo_main_loop_simulation = Simulation.from_config(
        montecarlo_main_loop_config,
        atom_data=atomic_dataset,
        virtual_packet_logging=False,
    )
    montecarlo_main_loop_simulation.run_convergence()
    montecarlo_main_loop_simulation.run_final()

    expected_hdf_store = regression_data.sync_hdf_store(
        montecarlo_main_loop_simulation
    )

    # Load compare data from refdata

    expected_nu = expected_hdf_store["/simulation/transport/output_nu"]
    expected_energy = expected_hdf_store["/simulation/transport/output_energy"]
    expected_nu_bar_estimator = expected_hdf_store[
        "/simulation/transport/nu_bar_estimator"
    ]
    expected_j_estimator = expected_hdf_store[
        "/simulation/transport/j_estimator"
    ]
    expected_hdf_store.close()
    actual_energy = montecarlo_main_loop_simulation.transport.output_energy
    actual_nu = montecarlo_main_loop_simulation.transport.output_nu
    actual_nu_bar_estimator = (
        montecarlo_main_loop_simulation.transport.nu_bar_estimator
    )
    actual_j_estimator = montecarlo_main_loop_simulation.transport.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-13)


@pytest.mark.xfail
def test_montecarlo_main_loop_vpacket_log(
    montecarlo_main_loop_config,
    regression_data,
    atomic_dataset,
):
    atomic_dataset = deepcopy(atomic_dataset)
    montecarlo_main_loop_config.montecarlo.no_of_virtual_packets = 5

    montecarlo_main_loop_simulation = Simulation.from_config(
        montecarlo_main_loop_config,
        atom_data=atomic_dataset,
        virtual_packet_logging=True,
    )
    montecarlo_main_loop_simulation.run_convergence()
    montecarlo_main_loop_simulation.run_final()

    expected_hdf_store = regression_data.sync_hdf_store(
        montecarlo_main_loop_simulation
    )

    expected_nu = expected_hdf_store["/simulation/transport/output_nu"]
    expected_energy = expected_hdf_store["/simulation/transport/output_energy"]
    expected_nu_bar_estimator = expected_hdf_store[
        "/simulation/transport/nu_bar_estimator"
    ]
    expected_j_estimator = expected_hdf_store[
        "/simulation/transport/j_estimator"
    ]
    expected_vpacket_log_nus = expected_hdf_store[
        "/simulation/transport/virt_packet_nus"
    ]
    expected_vpacket_log_energies = expected_hdf_store[
        "/simulation/transport/virt_packet_energies"
    ]

    actual_energy = (
        sim.transport.transport_state.packet_collection.output_energies
    )
    actual_nu = sim.transport.transport_state.packet_collection.output_nus
    actual_nu_bar_estimator = (
        sim.transport.transport_state.estimators.nu_bar_estimator
    )
    actual_j_estimator = sim.transport.transport_state.estimators.j_estimator
    expected_hdf_store.close()
    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy.value, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu.value, expected_nu, rtol=1e-13)
    npt.assert_allclose(
        actual_vpacket_log_nus, expected_vpacket_log_nus, rtol=1e-13
    )
    npt.assert_allclose(
        actual_vpacket_log_energies, expected_vpacket_log_energies, rtol=1e-13
    )
