import pytest
import numpy.testing as npt

from tardis.energy_input.gamma_ray_interactions import (
    compton_scatter,
    pair_creation,
    scatter_type,
)

from tardis.energy_input.GXPhoton import GXPhotonStatus


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
    initial_mu = basic_gamma_ray.direction.mu

    pair_creation(basic_gamma_ray)

    npt.assert_almost_equal(basic_gamma_ray.energy, 511.0)
    assert basic_gamma_ray.direction.mu != initial_mu


@pytest.mark.parametrize(
    ["compton_opacity", "photoabsorption_opacity", "total_opacity", "expected"],
    [
        (1, 0, 1, GXPhotonStatus.COMPTON_SCATTER),
        (0, 1, 1, GXPhotonStatus.PHOTOABSORPTION),
        (0, 0, 1, GXPhotonStatus.PAIR_CREATION),
    ],
)
def test_scatter_type(
    compton_opacity, photoabsorption_opacity, total_opacity, expected
):
    actual = scatter_type(
        compton_opacity, photoabsorption_opacity, total_opacity
    )
    assert actual == expected
