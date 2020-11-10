import pytest
import numpy as np

from tardis.energy_input.base import GammaRay, SphericalVector

@pytest.fixture(scope="function")
def basic_gamma_ray():
    return GammaRay(
        location = SphericalVector(1.e15, 0.5),
        direction = SphericalVector(1, 0.25),
        energy = 1000.e3,
        status = "InProcess",
        shell = 1
    )