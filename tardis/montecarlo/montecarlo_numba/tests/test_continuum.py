from copy import deepcopy

import numpy.testing as npt
import pytest

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
