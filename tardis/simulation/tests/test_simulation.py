import os

import pytest
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation

import pandas as pd
import pandas.util.testing as pdt


@pytest.fixture(scope='module')
def refdata(tardis_ref_data):
    def get_ref_data(key):
        return tardis_ref_data[os.path.join(
                'test_simulation', key)]
    return get_ref_data


@pytest.fixture(scope='module')
def config():
    return Configuration.from_yaml(
            'tardis/io/tests/data/tardis_configv1_verysimple.yml')


@pytest.fixture(scope='module')
def simulation_one_loop(
        atomic_data_fname, config,
        tardis_ref_data, generate_reference):
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    config.montecarlo.no_of_packets = int(4e4)
    config.montecarlo.last_no_of_packets = int(4e4)

    simulation = Simulation.from_config(config)
    simulation.run()

    if not generate_reference:
        return simulation
    else:
        simulation.hdf_properties = [
            'iterations_w',
            'iterations_t_rad',
            'iterations_electron_densities',
            'iterations_t_inner',
        ]
        simulation.model.hdf_properties = [
                't_radiative',
                'dilution_factor'
                ]
        simulation.runner.hdf_properties = [
                'j_estimator',
                'nu_bar_estimator',
                'output_nu',
                'output_energy'
                ]
        simulation.to_hdf(
                tardis_ref_data,
                '',
                'test_simulation'
        )
        simulation.model.to_hdf(
                tardis_ref_data,
                '',
                'test_simulation')
        simulation.runner.to_hdf(
                tardis_ref_data,
                '',
                'test_simulation')
        pytest.skip(
                'Reference data was generated during this run.')


@pytest.mark.parametrize('name', [
    'nu_bar_estimator', 'j_estimator', 't_radiative', 'dilution_factor',
    'output_nu', 'output_energy'
    ])
def test_plasma_estimates(
        simulation_one_loop, refdata, name):
    try:
        actual = getattr(
                simulation_one_loop.runner, name)
    except AttributeError:
        actual = getattr(
                simulation_one_loop.model, name)

    actual = pd.Series(actual)

    pdt.assert_almost_equal(
            actual,
            refdata(name)
            )


@pytest.mark.parametrize('name', [
    'iterations_w', 'iterations_t_rad',
    'iterations_electron_densities', 'iterations_t_inner'
    ])
def test_plasma_state_iterations(
        simulation_one_loop, refdata, name):
    actual = getattr(
        simulation_one_loop, name)

    try:
        actual = pd.Series(actual)
    except Exception:
        actual = pd.DataFrame(actual)

    pdt.assert_almost_equal(
            actual,
            refdata(name)
            )


#     assert_quantity_allclose(
#             t_rad, simulation_compare_data['test1/t_rad'] * u.Unit('K'), atol=0.0 * u.Unit('K'))
