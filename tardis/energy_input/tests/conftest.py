import pytest

from tardis.energy_input.util import SphericalVector, GXPhoton, GXPhotonStatus


@pytest.fixture(scope="function")
def basic_gamma_ray():
    """basic gamma ray fixture

    Returns
    -------
    GXPhoton
    """
    return GXPhoton(
        location=SphericalVector(1.0e15, 0.5),
        direction=SphericalVector(1, 0.25),
        energy=1000.0e3,
        status=GXPhotonStatus.IN_PROCESS,
        shell=1,
    )
