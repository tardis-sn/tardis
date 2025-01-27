import os

import astropy.units as u
import pytest

import tardis
from tardis.io.configuration.config_reader import Configuration
from tardis.simulation import Simulation

DATA_PATH = os.path.join(tardis.__path__[0], "plasma", "tests", "data")


@pytest.fixture
def config_init_trad_fname():
    return os.path.join(DATA_PATH, "config_init_trad.yml")


@pytest.mark.parametrize(
    "v_inner_boundary, v_outer_boundary",
    [
        (3350, 3650),
        (2900, 3750),
        (2900, 3850),
        (2900, 3900),
        (2950, 3750),
        (2950, 3850),
        (2950, 3900),
        (3050, 3750),
        (3050, 3850),
        (3050, 3900),
        (3150, 3750),
        (3150, 3850),
        (3150, 3900),
    ],
)
def test_plasma_vboundary(
    config_init_trad_fname,
    v_inner_boundary,
    v_outer_boundary,
    atomic_data_fname,
):
    tardis_config = Configuration.from_yaml(config_init_trad_fname)
    tardis_config.atom_data = atomic_data_fname
    tardis_config.model.structure.v_inner_boundary = (
        v_inner_boundary * u.km / u.s
    )
    tardis_config.model.structure.v_outer_boundary = (
        v_outer_boundary * u.km / u.s
    )
    simulation = Simulation.from_config(tardis_config)
