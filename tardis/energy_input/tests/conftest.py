import pytest

from tardis.energy_input.GXPhoton import GXPhoton, GXPhotonStatus


@pytest.fixture(scope="function")
def basic_gamma_ray():
    """basic gamma ray fixture

    Returns
    -------
    GXPhoton
    """
    return GXPhoton(
        location_r=1.0e15,
        location_theta=0.5,
        location_phi=0,
        direction_theta=1,
        direction_phi=0.25,
        energy=1000.0e3,
        status=GXPhotonStatus.IN_PROCESS,
        shell=1,
        activity=0,
    )
