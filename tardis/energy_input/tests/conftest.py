import pytest
import numpy as np

from tardis.energy_input.util import SphericalVector, GXPacket, GXPacketStatus


@pytest.fixture(scope="function")
def basic_gamma_ray():
    """basic gamma ray fixture

    Returns
    -------
    GXPacket
    """
    return GXPacket(
        location=SphericalVector(1.0e15, 0.5),
        direction=SphericalVector(1, 0.25),
        energy=1000.0e3,
        status=GXPacketStatus.IN_PROCESS,
        shell=1,
    )
