from copy import deepcopy

import numpy.testing as npt

from tardis.simulation import Simulation


def test_montecarlo_continuum(
    continuum_config,
    regression_data,
    nlte_atomic_dataset,
):
    nlte_atomic_dataset = deepcopy(nlte_atomic_dataset)
    continuum_simulation = Simulation.from_config(
        continuum_config,
        atom_data=nlte_atomic_dataset,
        virtual_packet_logging=False,
    )
    # continuum_simulation.run_convergence()

    continuum_simulation.run_final()

    expected_hdf_store = regression_data.sync_hdf_store(continuum_simulation)

    expected_nu = expected_hdf_store[
        "/simulation/transport/transport_state/output_nu"
    ]
    expected_energy = expected_hdf_store[
        "/simulation/transport/transport_state/output_energy"
    ]
    expected_hdf_store.close()

    transport_state = continuum_simulation.transport.transport_state
    actual_energy = transport_state.packet_collection.output_energies
    actual_nu = transport_state.packet_collection.output_nus

    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-13)


def test_run_iip_method(
    continuum_config,
    nlte_atomic_dataset,
):
    """Test that run_iip() method works directly."""
    nlte_atomic_dataset = deepcopy(nlte_atomic_dataset)
    continuum_simulation = Simulation.from_config(
        continuum_config,
        atom_data=nlte_atomic_dataset,
        virtual_packet_logging=False,
    )

    # Get transport solver and state
    transport_solver = continuum_simulation.transport
    transport_state = transport_solver.initialize_transport_state(
        continuum_simulation.simulation_state,
        continuum_simulation.plasma,
    )

    # Explicitly call run_iip() method
    v_packets_energy_hist = transport_solver.run_iip(
        transport_state,
        show_progress_bars=False,
    )

    # Verify IIP mode results
    assert v_packets_energy_hist is not None
    assert hasattr(transport_state, "estimators_continuum")
    assert transport_state.estimators_continuum is not None
    assert hasattr(transport_state.estimators_continuum, "j_blue_estimator")
    assert hasattr(
        transport_state.estimators_continuum, "alpha_stimulated_estimator"
    )
