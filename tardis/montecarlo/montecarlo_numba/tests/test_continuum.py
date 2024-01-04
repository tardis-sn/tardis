from copy import deepcopy

import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation


@pytest.fixture(scope="function")
def continuum_config(
    tardis_config_verysimple_nlte,
):
    continuum_config = Configuration.from_config_dict(
        tardis_config_verysimple_nlte
    )
    continuum_config.plasma.continuum_interaction.species = ["H I"]
    """"
    montecarlo_configuration.LEGACY_MODE_ENABLED = True
    # Setup model config from verysimple

    config_montecarlo_1e5_verysimple.montecarlo.last_no_of_packets = 1e5
    config_montecarlo_1e5_verysimple.montecarlo.no_of_virtual_packets = 0
    config_montecarlo_1e5_verysimple.montecarlo.iterations = 1
    config_montecarlo_1e5_verysimple.plasma.line_interaction_type = "macroatom"

    del config_montecarlo_1e5_verysimple["config_dirname"]
    """
    return continuum_config


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
    continuum_simulation.run_convergence()

    continuum_simulation.run_final()
    """
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
    transport_state = montecarlo_main_loop_simulation.transport.transport_state
    actual_energy = transport_state.packet_collection.output_energies
    actual_nu = transport_state.packet_collection.output_nus
    actual_nu_bar_estimator = transport_state.estimators.nu_bar_estimator
    actual_j_estimator = transport_state.estimators.j_estimator

    # Compare
    npt.assert_allclose(
        actual_nu_bar_estimator, expected_nu_bar_estimator, rtol=1e-13
    )
    npt.assert_allclose(actual_j_estimator, expected_j_estimator, rtol=1e-13)
    npt.assert_allclose(actual_energy, expected_energy, rtol=1e-13)
    npt.assert_allclose(actual_nu, expected_nu, rtol=1e-13)
"""
