from copy import deepcopy

import pytest

from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation

import numpy.testing as npt


@pytest.fixture(scope="function")
def continuum_config(
    tardis_config_verysimple_nlte,
):
    continuum_config = Configuration.from_config_dict(
        tardis_config_verysimple_nlte
    )
    continuum_config.plasma.continuum_interaction.species = ["H I"]
    continuum_config.plasma.nlte_ionization_species = []

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
