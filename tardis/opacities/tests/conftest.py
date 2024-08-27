import pytest
from copy import deepcopy
from tardis.simulation import Simulation


@pytest.fixture(scope="module")
def nb_simulation_verysimple(config_verysimple, atomic_dataset):
    atomic_data = deepcopy(atomic_dataset)
    sim = Simulation.from_config(config_verysimple, atom_data=atomic_data)
    sim.iterate(10)
    return sim
