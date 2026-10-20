from copy import deepcopy

import numpy.testing as npt

from tardis.simulation import Simulation


def test_montecarlo_transport_weighted(
    montecarlo_transport_config,
    montecarlo_transport_simulation,
    atomic_dataset,
    simple_weighted_packet_source,
):
    atomic_dataset = deepcopy(atomic_dataset)
    montecarlo_transport_simulation_weighted = Simulation.from_config(
        montecarlo_transport_config,
        atom_data=atomic_dataset,
        virtual_packet_logging=False,
        legacy_mode_enabled=True,
    )
    montecarlo_transport_simulation_weighted.packet_source = (
        simple_weighted_packet_source
    )
    montecarlo_transport_simulation_weighted.run_convergence()
    montecarlo_transport_simulation_weighted.run_final()

    expected_transport_state = (
        montecarlo_transport_simulation.transport.transport_state
    )
    transport_state = (
        montecarlo_transport_simulation_weighted.transport.transport_state
    )
    actual_energy = transport_state.packet_collection.output_energies
    actual_nu = transport_state.packet_collection.output_nus
    actual_nu_bar_estimator = transport_state.estimators_bulk.mean_frequency
    actual_j_estimator = transport_state.estimators_bulk.mean_intensity_total
    expected_energy = expected_transport_state.packet_collection.output_energies
    expected_nu = expected_transport_state.packet_collection.output_nus
    expected_nu_bar_estimator = (
        expected_transport_state.estimators_bulk.mean_frequency
    )
    expected_j_estimator = (
        expected_transport_state.estimators_bulk.mean_intensity_total
    )

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-14
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-14)
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-14)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-14)
