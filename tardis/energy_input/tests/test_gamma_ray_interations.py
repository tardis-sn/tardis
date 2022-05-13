import pytest
import numpy.testing as npt

from tardis.energy_input.gamma_ray_interactions import (
    compton_scatter,
    pair_creation_packet,
    scatter_type,
)
from tardis.energy_input.util import ELECTRON_MASS_ENERGY_KEV, H_CGS_KEV

from tardis.energy_input.GXPacket import GXPacketStatus


@pytest.mark.xfail(reason="To be implemented")
def test_get_compton_angle():
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_compton_scatter():
    assert False


def test_pair_creation(basic_gamma_ray):
    """

    Parameters
    ----------
    basic_gamma_ray : GammaRay object
    """
    initial_direction = basic_gamma_ray.direction

    pair_creation_packet(basic_gamma_ray)

    npt.assert_almost_equal(basic_gamma_ray.nu_cmf, ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV)
    assert basic_gamma_ray.direction != initial_direction


@pytest.mark.parametrize(
    ["compton_opacity", "photoabsorption_opacity", "total_opacity", "expected"],
    [
        (1, 0, 1, GXPacketStatus.COMPTON_SCATTER),
        (0, 1, 1, GXPacketStatus.PHOTOABSORPTION),
        (0, 0, 1, GXPacketStatus.PAIR_CREATION),
    ],
)
def test_scatter_type(
    compton_opacity, photoabsorption_opacity, total_opacity, expected
):
    actual = scatter_type(
        compton_opacity, photoabsorption_opacity, total_opacity
    )
    assert actual == expected
