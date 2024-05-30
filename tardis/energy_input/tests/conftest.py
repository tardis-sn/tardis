import numpy as np
import pytest

from tardis.energy_input.GXPacket import GXPacket, GXPacketStatus
from tardis.energy_input.util import H_CGS_KEV


@pytest.fixture(scope="function")
def basic_gamma_ray():
    """basic gamma ray fixture

    Returns
    -------
    GXPacket
    """
    return GXPacket(
        location=np.array([1.36375693e13, 4.10589818e14, 9.11718168e14]),
        direction=np.array([-0.97113853, 0.23134328, -0.05805379]),
        energy_rf=1e52,
        energy_cmf=1e52,
        nu_rf=1000.0e3 / H_CGS_KEV,
        nu_cmf=1000.0e3 / H_CGS_KEV,
        status=GXPacketStatus.IN_PROCESS,
        shell=1,
        time_current=1000,
    )
