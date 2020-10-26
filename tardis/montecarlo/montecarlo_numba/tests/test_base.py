import pytest
import pandas as pd
import os
import numpy.testing as npt
import numpy as np
import copy

import tardis.montecarlo.montecarlo_numba.base as base
from tardis.montecarlo.montecarlo_numba.numba_interface import PacketCollection
from tardis.montecarlo import montecarlo_configuration as montecarlo_configuration

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

@pytest.mark.xfail(reason='To be implemented')
def test_montecarlo_radial1d():
    assert False

#@pytest.mark.xfail(reason='To be implemented')
def test_montecarlo_main_loop(
    c_test_packet_collection, verysimple_numba_model, verysimple_numba_plasma, verysimple_estimators,
    nb_simulation_verysimple, set_seed_fixture
):
    montecarlo_configuration.single_packet_seed = 0
    
    output_packet_collection = PacketCollection(
        c_test_packet_collection.packets_input_nu,
        c_test_packet_collection.packets_input_mu,
        c_test_packet_collection.packets_input_energy,
        np.zeros(len(c_test_packet_collection.packets_input_nu), dtype=np.float64),
        np.zeros(len(c_test_packet_collection.packets_input_nu), dtype=np.float64)
    )

    set_seed_fixture(23111963)
    base.montecarlo_main_loop(
        output_packet_collection, verysimple_numba_model, verysimple_numba_plasma, verysimple_estimators,
        nb_simulation_verysimple.runner.spectrum_frequency.value, 0, [23111963]
    )

    npt.assert_allclose(output_packet_collection.packets_output_nu, c_test_packet_collection.packets_output_nu)
    npt.assert_allclose(output_packet_collection.packets_output_energy, c_test_packet_collection.packets_output_energy)