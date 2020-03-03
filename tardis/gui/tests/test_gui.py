import os

import pytest
from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import astropy.units as u


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

    returm simulation



def test_plasma_state_storer_store(simulation_one_loop)

    simulation = simulation_one_loop
