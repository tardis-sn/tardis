import pytest
import numpy.testing as npt

from tardis.energy_input.gamma_ray_interactions import compton_scatter, \
    pair_creation, photoabsorption, scatter_type

@pytest.mark.xfail(reason="To be implemented")
def test_compton_scatter():
    assert False


def test_pair_creation(basic_gamma_ray):
    initial_mu = basic_gamma_ray.direction.mu

    pair_creation(basic_gamma_ray)

    npt.assert_almost_equal(basic_gamma_ray.energy, 511.0)
    assert basic_gamma_ray.direction.mu != initial_mu


def test_photoabsorption(basic_gamma_ray):
    
    actual = photoabsorption(basic_gamma_ray)
    expected = basic_gamma_ray.energy

    npt.assert_almost_equal(actual, expected)
    assert basic_gamma_ray.status == "PhotoAbsorbed"


@pytest.mark.xfail(reason="To be implemented")
def test_scatter_type():
    assert False