from copy import deepcopy

import numpy.testing as npt
import pytest

from tardis.simulation import Simulation


@pytest.mark.xfail(reason="To be implemented")
def test_montecarlo_radial1d():
    raise AssertionError()


def test_montecarlo_transport(
    montecarlo_transport_config,
    regression_data,
    atomic_dataset,
):
    atomic_dataset = deepcopy(atomic_dataset)
    montecarlo_transport_simulation = Simulation.from_config(
        montecarlo_transport_config,
        atom_data=atomic_dataset,
        virtual_packet_logging=False,
        legacy_mode_enabled=True,
    )
    montecarlo_transport_simulation.run_convergence()
    montecarlo_transport_simulation.run_final()

    expected_hdf_store = regression_data.sync_hdf_store(
        montecarlo_transport_simulation
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
    transport_state = montecarlo_transport_simulation.transport.transport_state
    actual_energy = transport_state.packet_collection.output_energies
    actual_nu = transport_state.packet_collection.output_nus
    actual_nu_bar_estimator = transport_state.estimators_bulk.mean_frequency
    actual_j_estimator = transport_state.estimators_bulk.mean_intensity_total

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-13)


def test_montecarlo_transport_vpacket_log(
    montecarlo_transport_config,
    regression_data,
    atomic_dataset,
):
    atomic_dataset = deepcopy(atomic_dataset)
    montecarlo_transport_config.montecarlo.no_of_virtual_packets = 5

    montecarlo_transport_simulation = Simulation.from_config(
        montecarlo_transport_config,
        atom_data=atomic_dataset,
        virtual_packet_logging=True,
        legacy_mode_enabled=True,
    )
    montecarlo_transport_simulation.run_convergence()
    montecarlo_transport_simulation.run_final()

    transport = montecarlo_transport_simulation.transport

    assert transport.montecarlo_configuration.ENABLE_VPACKET_TRACKING is True

    expected_hdf_store = regression_data.sync_hdf_store(
        montecarlo_transport_simulation
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


def test_montecarlo_transport_hydrogen_only(
    config_verysimple_hydrogen_only, atomic_dataset, regression_data
):
    atomic_dataset = deepcopy(atomic_dataset)
    montecarlo_transport_simulation = Simulation.from_config(
        config_verysimple_hydrogen_only,
        atom_data=atomic_dataset,
        virtual_packet_logging=False,
    )
    montecarlo_transport_simulation.run_convergence()
    montecarlo_transport_simulation.run_final()

    transport = montecarlo_transport_simulation.transport.transport_state
    assert transport.j_estimator is not None
    expected_j_estimator = regression_data.sync_ndarray(transport.j_estimator)
    npt.assert_allclose(
        transport.j_estimator, expected_j_estimator, atol=0, rtol=1e-12
    )
