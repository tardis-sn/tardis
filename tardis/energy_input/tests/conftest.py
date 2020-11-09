import pytest
import numpy as np

from tardis.energy_input.base import GammaRay, SphericalVector

@pytest.fixture()
def basic_gamma_ray():
    return GammaRay(
        location = SphericalVector(1e15, np.pi/2),
        direction = SphericalVector(1, np.pi/4),
        energy = 1000.e3,
        status = "InProcess",
        shell = 1
    )