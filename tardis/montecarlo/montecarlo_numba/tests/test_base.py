import pytest
import pandas as pd
import os
import numpy.testing as npt
import numpy as np
import copy

import tardis.montecarlo.montecarlo_numba.base as base
from tardis.montecarlo.montecarlo_numba.numba_interface import PacketCollection, VPacketCollection
from tardis.montecarlo import montecarlo_configuration as montecarlo_configuration
import tardis.montecarlo.montecarlo_numba.numba_interface as numba_interface
import tardis.montecarlo.montecarlo_numba.numba_config as numba_config
import tardis.montecarlo.montecarlo_numba.r_packet as r_packet
import tardis.montecarlo.montecarlo_numba.single_packet_loop as spl

@pytest.fixture()
def read_c_test(tardis_ref_path):
    mode = "r"
    with pd.HDFStore(
        os.path.join(tardis_ref_path, "montecarlo_one_packet_compare_data.h5"), mode=mode
    ) as store:
        yield store

@pytest.fixture()
def c_test_packet_collection(read_c_test):
    input_data = read_c_test['/one_packet_loop']
    input_nu = input_data['input_nu'].values
    input_mu = input_data['input_mu'].values
    input_energy = input_data['input_energy'].values
    output_nu = input_data['output_nu'].values
    output_energy = input_data['output_energy'].values
    return PacketCollection(
        input_nu, input_mu, input_energy,
        output_nu, output_energy
    )

@pytest.fixture(scope="function")
def model():
    return numba_interface.NumbaModel(
        r_inner = np.array([6.912e14, 8.64e14], dtype=np.float64),
        r_outer = np.array([8.64e14, 1.0368e15], dtype=np.float64),
        time_explosion = 5.2e7
    )

@pytest.fixture(scope="function")
def estimators():
    return numba_interface.Estimators(
        j_estimator = np.array([0.0, 0.0], dtype=np.float64),
        nu_bar_estimator = np.array([0.0, 0.0], dtype=np.float64),
        j_blue_estimator = np.array([[1.e-10], [1.e-10]], dtype=np.float64),
        Edotlu_estimator = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    )

@pytest.fixture(scope="function")
def plasma():
    return numba_interface.NumbaPlasma(
        electron_density = np.array([1.0e9, 1.0e9], dtype=np.float64),
        line_list_nu = np.array([
            1.26318289e+16,
            1.26318289e+16,
            1.23357675e+16,
            1.23357675e+16,
            1.16961598e+16], dtype=np.float64),
        tau_sobolev = np.ones((500, 2), dtype=np.float64) * 1.0e-5,
        transition_probabilities = np.array([[0], [0]], dtype=np.float64),
        line2macro_level_upper = np.array([0, 0], dtype=np.int64),
        macro_block_references = np.array([0, 0], dtype=np.int64),
        transition_type = np.array([0, 0], dtype=np.int64),
        destination_level_id = np.array([0, 0], dtype=np.int64),
        transition_line_id = np.array([0, 0], dtype=np.int64)
    )


@pytest.mark.xfail(reason='To be implemented')
def test_montecarlo_radial1d():
    assert False

#@pytest.mark.xfail(reason='To be implemented')
def test_montecarlo_main_loop(
    c_test_packet_collection, model, plasma, estimators,
    nb_simulation_verysimple, set_seed_fixture, random_call_fixture
):
    montecarlo_configuration.single_packet_seed = 0
    
    output_packet_collection = PacketCollection(
        c_test_packet_collection.packets_input_nu,
        c_test_packet_collection.packets_input_mu,
        c_test_packet_collection.packets_input_energy,
        np.zeros(len(c_test_packet_collection.packets_input_nu), dtype=np.float64),
        np.zeros(len(c_test_packet_collection.packets_input_nu), dtype=np.float64)
    )

    numba_config.SIGMA_THOMSON = 6.652486e-25
    
    vpacket_collection = VPacketCollection(
        np.array([0, 0], dtype=np.float64), 0, np.inf, 
        0, 0
    )

    output_nus = np.empty_like(output_packet_collection.packets_output_nu)
    output_energies = np.empty_like(output_packet_collection.packets_output_nu)

    set_seed_fixture(23111963)
    seed = 23111963
    for i in range(len(c_test_packet_collection.packets_input_nu)):
        packet = r_packet.RPacket(model.r_inner[0],
                           output_packet_collection.packets_input_mu[i],
                           output_packet_collection.packets_input_nu[i],
                           output_packet_collection.packets_input_energy[i],
                           seed,
                           i, 0)

    
        spl.single_packet_loop(
            packet, model, plasma, estimators, vpacket_collection 
        )

        output_nus[i] = packet.nu

        if packet.status == r_packet.PacketStatus.REABSORBED:
            output_energies[i] = -packet.energy
        elif packet.status == r_packet.PacketStatus.EMITTED:
            output_energies[i] = packet.energy

        #random_call_fixture()

    # np.savetxt('scatter_output_energy.txt', output_energies)
    output_packet_collection.packets_output_energy[:] = output_energies[:]
    output_packet_collection.packets_output_nu[:] = output_nus[:]

    npt.assert_allclose(output_packet_collection.packets_output_nu, c_test_packet_collection.packets_output_nu, rtol=1e-12)
    npt.assert_allclose(output_packet_collection.packets_output_energy, c_test_packet_collection.packets_output_energy, rtol=1e-12)