import pytest
import numpy as np

from tardis.energy_input.base import main_gamma_ray_loop
from tardis.model import Radial1DModel
from tardis.io.config_reader import Configuration


@pytest.mark.skip(reason="To be implemented")
def test_initialize_photons():
    assert False


@pytest.mark.skip(reason="To be implemented")
def test_main_gamma_ray_loop():

    num_photons = 1e5
    num_photons = int(num_photons)

    np.random.seed(1)

    config = Configuration.from_yaml("tardis_configv1_verysimple.yml")

    model = Radial1DModel.from_config(config)

    ejecta_energy, ejecta_energy_r, escape_energy, radii = main_gamma_ray_loop(
        num_photons,
        model,
        path="decay_radiation.h5",
    )
