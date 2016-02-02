import numpy.testing as npt

import h5py
import pytest
from tardis.io import config_reader
from tardis.model import Radial1DModel
from tardis.simulation import Simulation


''' I don't get the purpose of this test. It looks like it's trying to do unit
tests for the Simulation module but what i see happening is a full tardis test
with only one iteration. As a result these tests will always break if something
in the Code changes.
To fix this, the test should be written in a way that their behavoir is
independent of other modules. If that cannot be achieved, move these tests to
test_tardis_full.py
'''


@pytest.fixture
def tardis_config(kurucz_atomic_data, tardis_config_verysimple):
    return config_reader.Configuration.from_config_dict(
        tardis_config_verysimple, atom_data=kurucz_atomic_data)


# For independent behavior the model should be read from a file, not generated
# at runtime
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

    npt.assert_allclose(
            t_rad, simulation_compare_data['test1/t_rad'], atol=0.0)
    npt.assert_allclose(w, simulation_compare_data['test1/w'], atol=0.0)


def test_packet_output(simulation_one_loop, simulation_compare_data):
    npt.assert_allclose(
            simulation_one_loop.runner.output_nu,
            simulation_compare_data['test1/output_nu'],
            atol=0.0)

    npt.assert_allclose(simulation_one_loop.runner.output_energy,
                        simulation_compare_data['test1/output_energy'],
                        atol=0.0)
