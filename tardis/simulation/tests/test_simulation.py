import os

import pytest
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import astropy.units as u


@pytest.fixture(scope="module")
def config():
    return Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
    )


@pytest.fixture(scope="module")
def simulation_one_loop(
    atomic_data_fname, config
):
    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    config.montecarlo.no_of_packets = int(4e4)
    config.montecarlo.last_no_of_packets = int(4e4)

    simulation = Simulation.from_config(config)
    simulation.run()
    return simulation

@pytest.mark.array_compare(file_format='pdhdf')
@pytest.mark.parametrize(
    "name",
    [
        "nu_bar_estimator",
        "j_estimator",
        "t_radiative",
        "dilution_factor",
        "output_nu",
        "output_energy",
    ],
)
def test_plasma_estimates(simulation_one_loop, name):
    try:
        actual = getattr(simulation_one_loop.runner, name)
    except AttributeError:
        actual = getattr(simulation_one_loop.model, name)


    actual = pd.Series(actual)
    return actual

@pytest.mark.array_compare(file_format='pdhdf')
@pytest.mark.parametrize(
    "name",
    [
        "iterations_w",
        "iterations_t_rad",
        "iterations_electron_densities",
        "iterations_t_inner",
    ],
)
def test_plasma_state_iterations(simulation_one_loop, name):
    actual = getattr(simulation_one_loop, name)

    try:
        actual = pd.Series(actual)
    except Exception:
        actual = pd.DataFrame(actual)

    return actual



@pytest.fixture(scope="module")
def simulation_without_loop(atomic_data_fname, config):

    config.atom_data = atomic_data_fname
    config.montecarlo.iterations = 2
    return Simulation.from_config(config)


def test_plasma_state_storer_store(
    atomic_data_fname, config, simulation_without_loop
):

    simulation = simulation_without_loop

    w_test = np.linspace(0, 1, 20)
    t_rad_test = np.linspace(12000, 9000, 20) * u.K
    electron_densities_test = pd.Series(np.linspace(1e7, 1e6, 20))
    t_inner_test = 12500 * u.K

    simulation.store_plasma_state(
        1, w_test, t_rad_test, electron_densities_test, t_inner_test
    )

    np.testing.assert_allclose(simulation.iterations_w[1, :], w_test)
    np.testing.assert_allclose(simulation.iterations_t_rad[1, :], t_rad_test)
    np.testing.assert_allclose(
        simulation.iterations_electron_densities[1, :], electron_densities_test
    )
    np.testing.assert_allclose(simulation.iterations_t_inner[1], t_inner_test)


def test_plasma_state_storer_reshape(
    atomic_data_fname, config, simulation_without_loop
):

    simulation = simulation_without_loop
    simulation.reshape_plasma_state_store(0)

    assert simulation.iterations_t_rad.shape == (1, 20)
    assert simulation.iterations_w.shape == (1, 20)
    assert simulation.iterations_electron_densities.shape == (1, 20)
    assert simulation.iterations_t_inner.shape == (1,)


#     assert_quantity_allclose(
#             t_rad, simulation_compare_data['test1/t_rad'] * u.Unit('K'), atol=0.0 * u.Unit('K'))
