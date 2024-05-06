from copy import deepcopy

import numpy.testing as npt
import pytest

from tardis.simulation import Simulation


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d():
    assert False


@pytest.fixture(scope="function")
def montecarlo_main_loop_config(
    config_montecarlo_1e5_verysimple,
):
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
        legacy_mode_enabled=True,
    )
    montecarlo_main_loop_simulation.run_convergence()
    montecarlo_main_loop_simulation.run_final()

    expected_hdf_store = regression_data.sync_hdf_store(
        montecarlo_main_loop_simulation
    )

    # Load compare data from refdata

    expected_nu = expected_hdf_store[
        "/simulation/transport/transport_state/output_nu"
    ]
    expected_energy = expected_hdf_store[
        "/simulation/transport/transport_state/output_energy"
    ]
    expected_nu_bar_estimator = expected_hdf_store[
        "/simulation/transport/transport_state/nu_bar_estimator"
    ]
    expected_j_estimator = expected_hdf_store[
        "/simulation/transport/transport_state/j_estimator"
    ]
    expected_hdf_store.close()
    transport_state = montecarlo_main_loop_simulation.transport.transport_state
    actual_energy = transport_state.packet_collection.output_energies
    actual_nu = transport_state.packet_collection.output_nus
    actual_nu_bar_estimator = (
        transport_state.radfield_mc_estimators.nu_bar_estimator
    )
    actual_j_estimator = transport_state.radfield_mc_estimators.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-13)


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
        legacy_mode_enabled=True,
    )
    montecarlo_main_loop_simulation.run_convergence()
    montecarlo_main_loop_simulation.run_final()

    transport = montecarlo_main_loop_simulation.transport

    assert transport.montecarlo_configuration.ENABLE_VPACKET_TRACKING is True

    expected_hdf_store = regression_data.sync_hdf_store(
        montecarlo_main_loop_simulation
    )

    expected_nu = expected_hdf_store[
        "/simulation/transport/transport_state/output_nu"
    ]
    expected_energy = expected_hdf_store[
        "/simulation/transport/transport_state/output_energy"
    ]
    expected_nu_bar_estimator = expected_hdf_store[
        "/simulation/transport/transport_state/nu_bar_estimator"
    ]
    expected_j_estimator = expected_hdf_store[
        "/simulation/transport/transport_state/j_estimator"
    ]
    expected_vpacket_log_nus = expected_hdf_store[
        "/simulation/transport/transport_state/virt_packet_nus"
    ]
    expected_vpacket_log_energies = expected_hdf_store[
        "/simulation/transport/transport_state/virt_packet_energies"
    ]

    transport_state = transport.transport_state

    actual_energy = transport_state.packet_collection.output_energies
    actual_nu = transport_state.packet_collection.output_nus
    actual_nu_bar_estimator = transport_state.nu_bar_estimator
    actual_j_estimator = transport_state.j_estimator
    actual_vpacket_log_nus = transport_state.vpacket_tracker.nus
    actual_vpacket_log_energies = transport_state.vpacket_tracker.energies

    expected_hdf_store.close()
    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator,
        expected_nu_bar_estimator,
        rtol=1e-12,
        atol=1e-15,
    )
    npt.assert_allclose(
        actual_j_estimator, expected_j_estimator, rtol=1e-12, atol=1e-15
    )
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-12, atol=1e-15)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-12, atol=1e-15)
    npt.assert_allclose(
        actual_vpacket_log_nus, expected_vpacket_log_nus, rtol=1e-12, atol=1e-15
    )
    npt.assert_allclose(
        actual_vpacket_log_energies,
        expected_vpacket_log_energies,
        rtol=1e-12,
        atol=1e-15,
    )
