import numpy.testing as npt

import numpy as np
import pandas as pd

import h5py
import pytest
from tardis.io import config_reader
from tardis.model import Radial1DModel
from tardis.simulation import Simulation
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

@pytest.fixture
def tardis_config(kurucz_atomic_data, tardis_config_verysimple):
    return config_reader.Configuration.from_config_dict(
        tardis_config_verysimple, atom_data=kurucz_atomic_data)


@pytest.fixture()
def raw_model(tardis_config):
    return Radial1DModel(tardis_config)


@pytest.fixture()
def simulation_one_loop(raw_model, tardis_config):
    sim = Simulation(tardis_config)
    sim.run_single_montecarlo(raw_model, 40000)

    return sim


@pytest.fixture()
def simulation_compare_data_fname():
    return 'tardis/simulation/tests/data/test_data.h5'

@pytest.fixture()
def simulation_mock_dataframes():
    lincol = ['atomic_number', 'ion_number', 'level_number_upper','wavelength']
    lin = np.hstack(( np.random.randint(0, 20, size = (5,3)),
                      np.random.random((5, 1)) ))
    lin = pd.DataFrame(lin,columns=lincol)
    lin.index.name = 'line_id'
    pass


@pytest.fixture()
def simulation_compare_data(simulation_compare_data_fname):
    return h5py.File(simulation_compare_data_fname, mode='r')


def test_plasma_estimates(simulation_one_loop, simulation_compare_data):
    t_rad, w = simulation_one_loop.runner.calculate_radiationfield_properties()

    npt.assert_allclose(simulation_one_loop.runner.nu_bar_estimator,
                        simulation_compare_data['test1/nubar_estimators'],
                        atol=0.0)
    npt.assert_allclose(simulation_one_loop.runner.j_estimator,
                        simulation_compare_data['test1/j_estimators'],
                        atol=0.0)

    assert_quantity_allclose(
            t_rad, simulation_compare_data['test1/t_rad'] * u.Unit('K'), atol=0.0 * u.Unit('K'))
    npt.assert_allclose(w, simulation_compare_data['test1/w'], atol=0.0)


def test_packet_output(simulation_one_loop, simulation_compare_data):
    assert_quantity_allclose(
            simulation_one_loop.runner.output_nu,
            simulation_compare_data['test1/output_nu'] * u.Unit('Hz'),
            atol=0.0 * u.Unit('Hz'))

    assert_quantity_allclose(simulation_one_loop.runner.output_energy,
                        simulation_compare_data['test1/output_energy'] * u.Unit('erg'),
                        atol=0.0 * u.Unit('erg'))


def test_make_source_function(simulation_one_loop,raw_model):
    # For now just test that the method can be called without exceptions
    tmp =  simulation_one_loop.make_source_function(raw_model)
