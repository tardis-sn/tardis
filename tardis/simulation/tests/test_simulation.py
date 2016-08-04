import copy

import h5py
import numpy.testing as npt
import pytest

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose, remote_data

from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
from tardis.simulation import Simulation


@pytest.fixture
def simulation_one_loop(atom_data, tardis_config_verysimple):
    tardis_config = Configuration.from_yaml(
        tardis_config_verysimple, atom_data=copy.deepcopy(atom_data)
    )
    raw_model = Radial1DModel(tardis_config)
    sim = Simulation(tardis_config)
    sim.run_single_montecarlo(raw_model, 40000)

    return sim


@pytest.fixture()
def simulation_compare_data_fname():
    return 'tardis/simulation/tests/data/test_data.h5'


@pytest.fixture()
def simulation_compare_data(simulation_compare_data_fname):
    return h5py.File(simulation_compare_data_fname, mode='r')


@remote_data
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


@remote_data
def test_packet_output(simulation_one_loop, simulation_compare_data):
    assert_quantity_allclose(
            simulation_one_loop.runner.output_nu,
            simulation_compare_data['test1/output_nu'] * u.Unit('Hz'),
            atol=0.0 * u.Unit('Hz'))

    assert_quantity_allclose(simulation_one_loop.runner.output_energy,
                        simulation_compare_data['test1/output_energy'] * u.Unit('erg'),
                        atol=0.0 * u.Unit('erg'))
