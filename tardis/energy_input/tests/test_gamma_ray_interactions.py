import pytest
import numpy as np
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
    """Test the Compton angle calculation"""
    assert False


@pytest.mark.xfail(reason="To be implemented")
def test_compton_scatter():
    """Test Compton scattering"""
    assert False


def test_pair_creation(basic_gamma_ray):
    """

    Parameters
    ----------
    basic_gamma_ray : GammaRay object
    """
    np.random.seed(2)

    initial_direction = basic_gamma_ray.direction
    basic_gamma_ray.nu_cmf = 2 * ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV

    pair_creation_packet(basic_gamma_ray)

    npt.assert_almost_equal(
        basic_gamma_ray.nu_cmf, ELECTRON_MASS_ENERGY_KEV / H_CGS_KEV
    )
    for i, j in zip(initial_direction, basic_gamma_ray.direction):
        assert i != j


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
    """Test the scattering type

    Parameters
    ----------
    compton_opacity : float
    photoabsorption_opacity : float
    total_opacity : float
    expected : list
        Expected parameters
    """
    actual = scatter_type(
        compton_opacity, photoabsorption_opacity, total_opacity
    )
    assert actual == expected
